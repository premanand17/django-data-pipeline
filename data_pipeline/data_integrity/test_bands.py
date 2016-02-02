''' Data integrity tests for cytobands index '''
from django.test import TestCase
from elastic.query import Query
from elastic.search import ElasticQuery, Search
from elastic.elastic_settings import ElasticSettings
import logging
from elastic.utils import ElasticUtils

logger = logging.getLogger(__name__)


class BandTest(TestCase):
    ''' Band tests. '''

    def test_chrom(self):
        ''' Check correct number of chromosomes. '''
        ids = ['X', 'Y']
        for i in range(22):
            ids.append(i+1)
        idx = ElasticSettings.idx('BAND', idx_type='CHROM')
        docs = Search(ElasticQuery(Query.ids(ids)), idx=idx, size=len(ids)).search().docs
        self.assertEqual(len(ids), len(docs), 'Check for chromosomes')
        for doc in docs:
            self.assertGreater(getattr(doc, 'length'), 1000000, 'Chromosome length')

    def test_bands(self):
        ''' Test cytobands type populated. '''
        (idx, idx_type) = ElasticSettings.idx('BAND', idx_type='BAND').split('/')
        self.assertGreater(ElasticUtils.get_docs_count(idx, idx_type), 1200)
