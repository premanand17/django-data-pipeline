import logging
from builtins import classmethod
from elastic.search import ElasticQuery, Search
from elastic.query import Query, TermsFilter
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
import json
logger = logging.getLogger(__name__)


class Gene(object):

    @classmethod
    def gene_mapping(cls, idx, idx_type):
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
             .add_property("nomenclature_authority_symbol", "string") \
             .add_property("dbxrefs", "object") \
             .add_property("pmids", "object") \
             .add_property("protein_id", "string")

        ''' create index and add mapping '''
        load = Loader()
        options = {"indexName": idx, "shards": 5}
        load.mapping(props, 'gene', analyzer=Loader.KEYWORD_ANALYZER, **options)

    @classmethod
    def ensembl_gene_parse(cls, ensembl_gene_parse):
        ''' Parse gene GTF file from ensembl to create gene index. '''
        gene_list = {}
        gene_list['docs'] = []
        for egene in ensembl_gene_parse:
            if egene.startswith('#'):
                continue
            parts = egene.split('\t')
            if 'ensembl' in parts[1] and parts[2] == 'gene':
                gi = {}
                gi['chromosome'] = parts[0]
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
    def gene2ensembl_parse(cls, gene2ens, idx):
        ''' Parse gene2ensembl file from NCBI and add entrez and protein_id's to gene index. '''
        genes = {}
        for gene in gene2ens:
            if gene.startswith('9606\t'):
                parts = gene.split('\t')
                gene_id = parts[1]
                ens_id = parts[2]
                prot_acc = parts[5]
                if ens_id in genes:
                    if prot_acc != '-':
                        if 'protein_id' in genes[ens_id]:
                            genes[ens_id]['protein_id'].append(prot_acc)
                        else:
                            genes[ens_id].update({'protein_id': [prot_acc]})
                else:
                    genes[ens_id] = {'dbxrefs': {'entrez': gene_id}}
                    if prot_acc != '-':
                        genes[ens_id].update({'protein_id': [prot_acc]})

        query = ElasticQuery(Query.ids(list(genes.keys())))
        docs = Search(query, idx=idx, size=80000).search().docs

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
                genes[parts[1]]["PMID"].append(pmid)
            else:
                genes[parts[1]] = {"PMID": [pmid]}
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
            db = dbx[0].replace('HGNC:HGNC', 'HGNC')
            arr[db.lower()] = dbx[1]
        gi['dbxrefs'] = arr
