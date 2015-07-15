''' Used to generate the gene documents for indexing. '''

import logging
from builtins import classmethod
import sys
import csv
import re
import zipfile
import os
from elastic.search import ElasticQuery, Search
from elastic.query import Query, TermsFilter
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
import json

logger = logging.getLogger(__name__)


class Gene(object):
    ''' Gene class to define functions for building gene related indices. '''

    @classmethod
    def gene_mapping(cls, idx, idx_type):
        ''' Load the mapping for the gene index. '''
        props = MappingProperties(idx_type)
        props.add_property("symbol", "string", analyzer="full_name") \
             .add_property("synonyms", "string", analyzer="full_name") \
             .add_property("chromosome", "string") \
             .add_property("source", "string") \
             .add_property("start", "integer") \
             .add_property("end", "integer") \
             .add_property("strand", "string") \
             .add_property("description", "string") \
             .add_property("biotype", "string") \
             .add_property("dbxrefs", "object") \
             .add_property("pmids", "string")

        ''' create index and add mapping '''
        load = Loader()
        options = {"indexName": idx, "shards": 5}
        load.mapping(props, 'gene', analyzer=Loader.KEYWORD_ANALYZER, **options)

    @classmethod
    def ensembl_gene_parse(cls, ensembl_gene_parse):
        ''' Parse gene GTF file from ensembl to create gene index. '''
        gene_list = {}
        gene_list['docs'] = []
        allowed_chr = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                       "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                       "21", "22", "X", "Y"]
        for egene in ensembl_gene_parse:
            if egene.startswith('#'):
                continue
            parts = egene.split('\t')
            if parts[2] == 'gene':
                gi = {}
                gi['chromosome'] = parts[0].upper()
                if gi['chromosome'] not in allowed_chr:
                    continue
                gi['source'] = parts[1]
                gi['start'] = int(parts[3])
                gi['stop'] = int(parts[4])
                gi['strand'] = parts[6]

                attrs = parts[8].split(';')
                for attr in attrs:
                    a = attr.strip().split(" ")
                    if a[0] == 'gene_id':
                        gi['_id'] = a[1][1:-1]
                    elif a[0] == 'gene_biotype':
                        gi['biotype'] = a[1][1:-1]
                    elif a[0] == 'gene_name':
                        gi['symbol'] = a[1][1:-1]
                if '_id' in gi:
                    gene_list['docs'].append(gi)
        return gene_list

    @classmethod
    def gene2ensembl_parse(cls, gene2ens, idx, idx_type):
        ''' Parse gene2ensembl file from NCBI and add entrez to gene index. '''
        genes = {}
        for gene in gene2ens:
            if gene.startswith('9606\t'):
                parts = gene.split('\t')
                gene_id = parts[1]
                ens_id = parts[2]
#                 prot_acc = parts[5]
                if ens_id not in genes:
                    genes[ens_id] = {'dbxrefs': {'entrez': gene_id}}

        query = ElasticQuery(Query.ids(list(genes.keys())))
        docs = Search(query, idx=idx, idx_type=idx_type, size=80000).search().docs

        chunk_size = 450
        for i in range(0, len(docs), chunk_size):
            docs_chunk = docs[i:i+chunk_size]
            json_data = ''
            for doc in docs_chunk:
                ens_id = doc._meta['_id']
                idx_type = doc.type()
                doc_data = {"update": {"_id": ens_id, "_type": idx_type,
                                       "_index": idx, "_retry_on_conflict": 3}}
                json_data += json.dumps(doc_data) + '\n'
                json_data += json.dumps({'doc': genes[ens_id]}) + '\n'
            if json_data != '':
                Loader().bulk_load(idx, idx_type, json_data)

    @classmethod
    def ensmart_gene_parse(cls, ensmart_f, idx, idx_type):
        ''' For those gene docs missing a dbxrefs.entrez use Ensembl Mart to
        fill in. '''
        genes = {}
        for ensmart in ensmart_f:
            parts = ensmart.split('\t')
            ens_id = parts[0]
            gene_id = parts[1]
            swissprot = parts[2].strip()
            trembl = parts[3].strip()
            if gene_id == '':
                continue
            if ens_id in genes:
                if genes[ens_id]['dbxrefs']['entrez'] != gene_id:
                    genes[ens_id]['dbxrefs']['entrez'] = None
                else:
                    if swissprot != '':
                        cls._add_to_dbxref(genes[ens_id], 'swissprot', swissprot)
                    if trembl != '':
                        cls._add_to_dbxref(genes[ens_id], 'trembl', trembl)
            else:
                genes[ens_id] = {'dbxrefs': {'entrez': gene_id}}
                if swissprot != '':
                    genes[ens_id]['dbxrefs'].update({'swissprot': swissprot})
                if trembl != '':
                    genes[ens_id]['dbxrefs'].update({'trembl': trembl})

        '''  search for the entrez ids '''
        query = ElasticQuery(Query.ids(list(genes.keys())))
        docs = Search(query, idx=idx, idx_type=idx_type, size=80000).search().docs
        chunk_size = 450
        for i in range(0, len(docs), chunk_size):
            docs_chunk = docs[i:i+chunk_size]
            json_data = ''
            for doc in docs_chunk:
                ens_id = doc._meta['_id']
                if 'dbxrefs' in doc.__dict__:
                    dbxrefs = getattr(doc, 'dbxrefs')
                else:
                    dbxrefs = {}

                if ('entrez' in genes[ens_id]['dbxrefs'] and
                    'entrez' in dbxrefs and
                   dbxrefs['entrez'] != genes[ens_id]['dbxrefs']['entrez']):
                    logger.warn('Multiple entrez ids for ensembl id: '+ens_id)
                    continue

                idx_type = doc.type()
                doc_data = {"update": {"_id": ens_id, "_type": idx_type,
                                       "_index": idx, "_retry_on_conflict": 3}}
                json_data += json.dumps(doc_data) + '\n'
                json_data += json.dumps({'doc': genes[ens_id]}) + '\n'
            if json_data != '':
                Loader().bulk_load(idx, idx_type, json_data)

    @classmethod
    def _add_to_dbxref(cls, gene, db, dbxref):
        if db in gene['dbxrefs']:
            if not isinstance(gene['dbxrefs'][db], list):
                if dbxref == gene['dbxrefs'][db]:
                    return
                gene['dbxrefs'][db] = [gene['dbxrefs'][db]]
            elif dbxref in gene['dbxrefs'][db]:
                return
            gene['dbxrefs'][db].append(dbxref)
        else:
            gene['dbxrefs'].update({db: dbxref})

    @classmethod
    def gene_info_parse(cls, gene_infos, idx):
        ''' Parse gene_info file from NCBI and add info to gene index. '''

        # tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description type_of_gene
        # Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority Nomenclature_status
        # Other_designations Modification_date]
        genes = {}
        for gene_info in gene_infos:
            if gene_info.startswith('9606\t'):
                parts = gene_info.split('\t')
                gene = {"synonyms": parts[4].split("|")}
                cls._set_dbxrefs(parts[1], parts[5], gene)
                gene.update({"description": parts[8]})
                genes[parts[1]] = gene
        cls._update_gene(genes, idx)

    @classmethod
    def gene_pub_parse(cls, gene_pubs, idx):
        ''' Parse gene2pubmed file from NCBI and add PMIDs to gene index. '''
        genes = {}
        for gene_pub in gene_pubs:
            if not gene_pub.startswith('9606\t'):
                continue
            parts = gene_pub.split('\t')
            pmid = parts[2].strip()
            if parts[1] in genes:
                genes[parts[1]]["pmids"].append(pmid)
            else:
                genes[parts[1]] = {"pmids": [pmid]}
        cls._update_gene(genes, idx)

    @classmethod
    def _update_gene(cls, genes, idx):
        ''' Use genes data to update the index. '''
        gene_keys = list(genes.keys())
        chunk_size = 450
        for i in range(0, len(genes), chunk_size):
            chunk_gene_keys = gene_keys[i:i+chunk_size]
            json_data = ''

            query = ElasticQuery.filtered(Query.match_all(),
                                          TermsFilter.get_terms_filter("dbxrefs.entrez", chunk_gene_keys))
            docs = Search(query, idx=idx, size=chunk_size).search().docs
            for doc in docs:
                ens_id = doc._meta['_id']
                idx_type = doc.type()
                entrez = getattr(doc, 'dbxrefs')['entrez']
                doc_data = {"update": {"_id": ens_id, "_type": idx_type,
                                       "_index": idx, "_retry_on_conflict": 3}}
                json_data += json.dumps(doc_data) + '\n'
                json_data += json.dumps({'doc': genes[entrez]}) + '\n'
            if json_data != '':
                Loader().bulk_load(idx, idx_type, json_data)

    @classmethod
    def _set_dbxrefs(cls, entrez, dbxrefs, gi):
        arr = {'entrez': entrez}
        if dbxrefs == '-':
            return arr

        for dbxref in dbxrefs.split('|'):
            dbx = dbxref.rsplit(":", 1)
            if len(dbx) != 2:
                logger.warn('DBXREF PARSE: '+dbxref)
                continue
            dbx[0] = dbx[0].lower()
            if dbx[0] == 'ensembl':
                continue
            db = dbx[0].replace('hgnc:hgnc', 'hgnc')
            arr[db] = dbx[1]
        gi['dbxrefs'] = arr

    @classmethod
    def gene_interaction_parse(cls, download_file, stage_output_file, section):
        cls._psimitab(download_file, stage_output_file, section)

    @classmethod
    def _psimitab(cls, download_file, stage_output_file, section):
        abs_path_download_dir = os.path.dirname(download_file)
        zf = zipfile.ZipFile(download_file, 'r')
        # print(zf.namelist())

        import_file_exists = False
        # target_path = '/dunwich/scratch/prem/tmp/download/intact.txt'

        if import_file_exists is not True:
            stage_output_file_handler = open(stage_output_file, 'w')
            header_line = 'interactorA' + '\t' + 'interactorB' + '\t' + 'pubmed' + '\n'

            stage_output_file_handler.write(header_line)

            if 'intact.txt' in zf.namelist():
                # base_dir_path = args[3]
                print('Extracting the zip file...')
                target_path = zf.extract(member='intact.txt', path=abs_path_download_dir)
                line_number = 0
                with open(target_path, encoding='utf-8') as csvfile:
                    reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
                    for row in reader:
                            # print(row)
                            # check for taxid
                            re.compile('taxid:9606')
                            is_human_A = cls._check_tax_id(row['Taxid interactor A'], 'taxid:9606')
                            is_human_B = cls._check_tax_id(row['Taxid interactor B'], 'taxid:9606')

                            if is_human_A and is_human_B:
                                pass
                            else:
                                continue

                            # check for pubmed id/evidence
                            cleaned_pubmed_id = cls._clean_id(row['Publication Identifier(s)'], 'pubmed:\d+')

                            if cleaned_pubmed_id is None:
                                cleaned_pubmed_id = ''
                            # xref id
                            cleaned_xref_id_A = cls._clean_id(row['Xref(s) interactor A'], 'ensembl:ENSG\d+')
                            cleaned_xref_id_B = cls._clean_id(row['Xref(s) interactor B'], 'ensembl:ENSG\d+')

                            if (cleaned_xref_id_A is not None and
                               cleaned_xref_id_B is not None and
                               cleaned_xref_id_A != cleaned_xref_id_B):
                                line_number += 1
                                line = cleaned_xref_id_A + '\t' + cleaned_xref_id_B + '\t' + cleaned_pubmed_id + '\n'
                                stage_output_file_handler.write(line)

                stage_output_file_handler.close()
                cls._process_interaction_out_file(stage_output_file, section)
        else:
            cls._process_interaction_out_file(stage_output_file, section)

    @classmethod
    def _process_interaction_out_file(cls, target_path, section):
        print('_process_interaction_out_file file called...start processing')

        dict_container = dict()
        line_number = 0

        with open(target_path) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader:
                line_number += 1
                sys.stdout.write('.')
                # print(row)
                cleaned_xref_id_A = row['interactorA']
                cleaned_xref_id_B = row['interactorB']
                cleaned_pubmed_id = row['pubmed']

                if cleaned_xref_id_A == cleaned_xref_id_B:
                    continue

                dict_interactorA = {'interactor': cleaned_xref_id_A, 'pubmed': cleaned_pubmed_id}
                dict_interactorB = {'interactor': cleaned_xref_id_B, 'pubmed': cleaned_pubmed_id}

                if cleaned_xref_id_A in dict_container:
                    current_A = dict_container[cleaned_xref_id_A]
                    if dict_interactorB not in current_A:
                        current_A.append(dict_interactorB)
                        dict_container[cleaned_xref_id_A] = current_A
                else:
                    dict_container[cleaned_xref_id_A] = [dict_interactorB]

                if cleaned_xref_id_B in dict_container:
                    current_B = dict_container[cleaned_xref_id_B]
                    if dict_interactorA not in current_B:
                        current_B.append(dict_interactorA)
                        dict_container[cleaned_xref_id_B] = current_B
                else:
                    dict_container[cleaned_xref_id_B] = [dict_interactorA]

            cls._create_json_output(dict_container, target_path, section)
            print('GENE INTERACTION STAGE COMPLETE')

    @classmethod
    def _create_json_output(cls, dict_container, target_file_path, section):
        count = 0
        dict_keys = dict_container.keys()
        print('No of docs to load ' + str(len(dict_keys)))
        json_target_file_path = target_file_path.replace(".out", ".json")
        print(json_target_file_path)
        # json_target_file_path = target_file_path + '.json'

        load_mapping = True
        with open(json_target_file_path, mode='w', encoding='utf-8') as f:
            f.write('{"docs":[\n')

            for gene in dict_container:
                int_object = dict()
                int_object["_id"] = gene
                int_object["_parent"] = gene
                int_object["interactors"] = dict_container[gene]
                int_object["interaction_source"] = section['source']
                f.write(json.dumps(int_object))
                count += 1

                if len(dict_keys) == count:
                    f.write('\n')
                else:
                    f.write(',\n')

            f.write('\n]}')
        logger.debug("No. genes to load "+str(count))
        logger.debug("Json written to " + json_target_file_path)
        logger.debug("Load mappings")

        if load_mapping:
            status = cls._load_interaction_mappings(section)
            print(str(status))

    @classmethod
    def _load_interaction_mappings(cls, section):
        interaction_mapping = MappingProperties("interactions", "gene")
        interaction_mapping.add_property("interactors", "object", index="not_analyzed")
        interaction_mapping.add_property("interaction_source", "string")
        load = Loader()
        idx = section['index']
        options = {"indexName": idx, "shards": 1}
        status = load.mapping(interaction_mapping, "interactions", **options)
        return status

    @classmethod
    def _check_tax_id(cls, search_str, idpattern):
        p = re.compile(idpattern)
        m = p.search(search_str)

        if m:
            return True
        else:
            return False

    @classmethod
    def _clean_id(cls, search_str, idpattern):
        p = re.compile(idpattern)
        m = p.search(search_str)
        matched_id = ''
        if m:
            # print('Match found: ', m.group())
            matched_id = m.group()
            split_id = matched_id.split(sep=':')
            if split_id:
                return split_id[1]
        else:
            pass

        return None
