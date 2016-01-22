''' Data integrity tests for region index '''
from django.test import TestCase
from elastic.elastic_settings import ElasticSettings
import logging
from elastic.query import Query
from elastic.search import Search, ElasticQuery
from data_pipeline.data_integrity.utils import DataIntegrityUtils
from region import utils
logger = logging.getLogger(__name__)


class RegionDataTest(TestCase):
    '''Region test '''
    IDX_KEY = 'REGION'
    IDX_TYPE_KEYS = ['STUDY_HITS', 'DISEASE_LOCUS', 'REGION']
    DOC_COUNTS = {
        'STUDY_HITS': 1760,
        'DISEASE_LOCUS': 900,
        'REGION': 342
    }

    def test_hit_attributes(self):
        '''Fetch random genes from elastic and compare the same with the results fetched via ensembl restful query'''

        for idx_type_key in RegionDataTest.IDX_TYPE_KEYS:
            idx = ElasticSettings.idx(RegionDataTest.IDX_KEY, idx_type_key)
            (idx, idx_type) = idx.split('/')

            docs = DataIntegrityUtils.get_rdm_docs(idx, idx_type, qbool=Query.match_all(), sources=[], size=1)

    def test_docs_count(self):
        '''Check the number of docs in a given index/index-type'''

        for idx_type_key in RegionDataTest.IDX_TYPE_KEYS:
            idx = ElasticSettings.idx(RegionDataTest.IDX_KEY, idx_type_key)
            (idx, idx_type) = idx.split('/')

            doc_count = DataIntegrityUtils.get_docs_count(idx, idx_type)
            self.assertEqual(doc_count, RegionDataTest.DOC_COUNTS[idx_type_key],
                             "Count of docs in the "+idx_type_key+" index are correct")

    def test_region_attributes(self):
        ''' test region attributes '''
        idx = ElasticSettings.idx(RegionDataTest.IDX_KEY, 'REGION')
        (idx, idx_type) = idx.split('/')
        docs = DataIntegrityUtils.get_rdm_docs(idx, idx_type, qbool=Query.match_all(), sources=[], size=1)
        newRegion = utils.Region.pad_region_doc(docs[0])

        if len(getattr(newRegion, "genes")) > 0:
            query = ElasticQuery(Query.ids(getattr(newRegion, "genes")))
            resultObject = Search(query, idx=ElasticSettings.idx('GENE', 'GENE'),
                                  size=len(getattr(newRegion, "genes"))).search()
            self.assertEqual(len(getattr(newRegion, "genes")), resultObject.hits_total,
                             "All genes on region found in GENE index")

        if len(getattr(newRegion, "studies")) > 0:
            query = ElasticQuery(Query.ids(getattr(newRegion, "studies")))
            resultObject = Search(query, idx=ElasticSettings.idx('STUDY', 'STUDY'),
                                  size=len(getattr(newRegion, "studies"))).search()
            self.assertEqual(len(getattr(newRegion, "studies")), resultObject.hits_total,
                             "All study ids for region found in STUDY index")

        if len(getattr(newRegion, "pmids")) > 0:
            query = ElasticQuery(Query.ids(getattr(newRegion, "pmids")))
            resultObject = Search(query, idx=ElasticSettings.idx('PUBLICATION', 'PUBLICATION'),
                                  size=len(getattr(newRegion, "pmids"))).search()
            self.assertEqual(len(getattr(newRegion, "pmids")), resultObject.hits_total,
                             "All PMIDs for region found in PUBLICATION index")
