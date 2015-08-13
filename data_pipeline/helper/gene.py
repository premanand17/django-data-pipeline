''' Used to generate the gene documents for indexing. '''

import logging
from builtins import classmethod
from elastic.search import ElasticQuery, Search
from elastic.query import Query, TermsFilter
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
import json

logger = logging.getLogger(__name__)


class Gene(object):
    ''' Gene class to define functions for building gene related indices.

    The gene index is built by parsing the following:
    1. Ensembl gene GTF (ensembl_id, symbol, biotype, chromosome, source, start, stop, strand)
    2. NCBI gene2ensembl (entrez id)
    3. Ensembl Mart (missing entrez ids, swissprot, trembl)
    4. NCBI gene_info (synonyms, dbxrefs, description)
    5. NCBI gene2pubmed (pmids)
    '''

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
             .add_property("pmids", "string") \
             .add_property("suggest", "completion",
                           index_analyzer="full_name", search_analyzer="full_name")

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
                suggests = parts[4].split("|")
                if 'dbxrefs' in gene:
                    suggests.extend(list(gene['dbxrefs'].values()))
                suggests.append(parts[2])
                gene['suggest'] = {}
                gene['suggest']["input"] = suggests
                gene['suggest']["weight"] = 50
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

            db = dbx[0].replace('hgnc:hgnc', 'hgnc')
            arr[db] = dbx[1]
        gi['dbxrefs'] = arr

    @classmethod
    def _convert_entrezid2ensembl(cls, gene_sets, section, log_output_file_handler=None, log_conversion=True):
        '''Converts given set of entrez ids to ensembl ids by querying the gene index dbxrefs'''
        query = ElasticQuery.filtered(Query.match_all(),
                                      TermsFilter.get_terms_filter("dbxrefs.entrez", gene_sets))
        docs = Search(query, idx=section['index'], size=1000000).search().docs
        ensembl_ids = []
        for doc in docs:
            ens_id = doc._meta['_id']
            ensembl_ids.append(ens_id)

        if log_conversion:
            if log_output_file_handler is not None:
                cls._log_entrezid2ensembl_coversion(gene_sets, ensembl_ids, log_output_file_handler)

        return ensembl_ids

    @classmethod
    def _log_entrezid2ensembl_coversion(cls, entrez_genes_in, ensembl_genes_out,  log_output_file_handler):
        '''Logs the conversion rates between entrez2ensembl and also stores the input entrez ids for later reference'''
        entrez_genes_count = len(entrez_genes_in)
        ensembl_genes_count = len(ensembl_genes_out)
        try:
            less_more = entrez_genes_count/ensembl_genes_count
        except ZeroDivisionError:
            less_more = 2

        if less_more > 1:
            diff = entrez_genes_count - ensembl_genes_count
            diff_text = "Less("+str(diff) + ")"
        elif less_more < 1:
            diff = ensembl_genes_count - entrez_genes_count
            diff_text = "More("+str(diff) + ")"
        elif less_more == 1:
            diff = entrez_genes_count - ensembl_genes_count
            diff_text = "Equal("+str(diff) + ")"

        log_output_file_handler.write(','.join(entrez_genes_in) + "\t" +
                                      ','.join(ensembl_genes_out) + "\t" + str(entrez_genes_count) + '/' +
                                      str(ensembl_genes_count) + "\t" + diff_text + "\n")
