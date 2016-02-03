''' Data integrity tests for hapmap recombination rates index '''
import logging

from django.test import TestCase

from elastic.elastic_settings import ElasticSettings
from elastic.utils import ElasticUtils


logger = logging.getLogger(__name__)


class RecombinationTest(TestCase):
    ''' HapMap Recombination tests. '''

    def test_data_loaded(self):
        ''' Test cytobands type populated. '''
        (idx, idx_type) = ElasticSettings.idx('HAPMAP', idx_type='HAPMAP').split('/')
        self.assertGreater(ElasticUtils.get_docs_count(idx, idx_type), 3000000)
