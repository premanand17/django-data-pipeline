''' Data integrity tests for gene index '''
from django.test import TestCase
from elastic.elastic_settings import ElasticSettings
import logging
from elastic.query import Query
from elastic.search import ElasticQuery, Search
import requests
import sys
import json
import re
from data_pipeline.data_integrity.utils import DataIntegrityUtils
logger = logging.getLogger(__name__)


class GeneDataTest(TestCase):
    '''Gene test '''

    def test_gene_attributes(self):
        '''Fetch random genes from elastic and compare the same with the results fetched via ensembl restful query'''
        idx_key = 'GENE'
        idx_type_key = 'GENE'

        idx = ElasticSettings.idx(idx_key, idx_type_key)
        (idx, idx_type) = idx.split('/')

        geneids = DataIntegrityUtils.get_rdm_feature_ids(idx, idx_type, qbool=Query.match_all(), sources=[],
                                                         field=None, size=2)

        logger.debug("===Testing genes===")
        logger.debug(geneids)
        query = ElasticQuery(Query.ids(geneids))
        elastic = Search(query, idx=ElasticSettings.idx(idx_key, idx_type_key), size=len(geneids))

        docs_by_geneid = elastic.search().docs

        # "_source":{"symbol": "RP11-376M2.2", "start": 42975689, "biotype": "sense_intronic", "chromosome": "17",
        # "source": "havana", "strand": "-", "stop": 42977275}
        for doc in docs_by_geneid:
            gene_id_pipeline = doc.doc_id()
            index_pipeline = doc.index()
            start_pipeline = getattr(doc, "start")
            stop_pipeline = getattr(doc, "stop")
            chromosome_pipeline = getattr(doc, "chromosome")

            biotype_pipeline = getattr(doc, "biotype")
            strand_pipeline = getattr(doc, "strand")
            strand_pipeline = -1 if strand_pipeline == '-' else 1
            symbol_pipeline = getattr(doc, "symbol")
            source_pipeline = getattr(doc, "source")

            # genes_hg38_v0.0.2
            pattern = re.compile('genes_\w\w(\d+)', re.IGNORECASE)
            match = pattern.match(index_pipeline)
            assembly_number_pipeline = None
            if match:
                assembly_number_pipeline = match.group(1)

            ensembl_gene_data = self.check_at_ensembl(gene_id_pipeline)

            pattern = re.compile('GRCh(\d+)', re.IGNORECASE)
            match = pattern.match(ensembl_gene_data['assembly_name'])

            assembly_number_ens = None
            if match:
                assembly_number_ens = match.group(1)

            self.assertEqual(assembly_number_pipeline, assembly_number_ens, "Assembly number is ok")
            self.assertEqual(gene_id_pipeline, ensembl_gene_data['id'], "Gene Id number is ok")
            self.assertEqual(start_pipeline, ensembl_gene_data['start'], "start is ok")
            self.assertEqual(stop_pipeline, ensembl_gene_data['end'], "stop is ok")
            self.assertEqual(chromosome_pipeline, ensembl_gene_data['seq_region_name'], "chr is ok")
            self.assertEqual(strand_pipeline, ensembl_gene_data['strand'], "strand is ok")

            self.assertEqual(biotype_pipeline, ensembl_gene_data['biotype'], "biotype is ok")
            self.assertEqual(symbol_pipeline, ensembl_gene_data['display_name'], "symbol/display_name is ok")
            self.assertEqual(source_pipeline, ensembl_gene_data['source'], "source is ok")

    def check_at_ensembl(self, gene_id):
        '''Lookup ensembl via restful call'''
        server = "http://rest.ensembl.org"
        ext = "/lookup/id/" + gene_id + "?content-type=application/json;expand=1;db_type=core;object_type=Gene"

        url = server+ext
        logger.debug(url)

        response = requests.get(url)

        if not response.ok:
            response.raise_for_status()
            sys.exit()

        data = json.loads(response.content.decode('utf-8'))
        return data
