''' Tests for the staging and loading. '''
from django.test import TestCase
from django.core.management import call_command
import os
import data_pipeline
import shutil
from data_pipeline.utils import IniParser
from elastic.elastic_settings import ElasticSettings
import requests
from elastic.search import Search
from data_pipeline.helper.gene import Gene

INI_FILE = os.path.join(os.path.dirname(__file__), 'test_download.ini')
TEST_DATA_DIR = os.path.dirname(data_pipeline.__file__) + '/tests/data'
INI_CONFIG = IniParser().read_ini(INI_FILE)


def tearDownModule():
    shutil.rmtree(TEST_DATA_DIR + '/STAGE')
    # remove index created
    requests.delete(ElasticSettings.url() + '/' + INI_CONFIG['GENE_HISTORY']['index'])


class StageTest(TestCase):

    def test_gene_interaction_parse(self):
        ''' Test ini file downloads. '''
        call_command('pipeline', '--steps', 'stage', sections='INTACT',
                     dir=TEST_DATA_DIR, ini=INI_FILE)

        self.assertTrue(os.path.exists(TEST_DATA_DIR + '/STAGE/INTACT'))
        self.assertTrue(os.path.isfile(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.out'))
        self.assertTrue(os.path.isfile(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.json'))
        self.assertTrue(os.access(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.json', os.R_OK))


class LoadTest(TestCase):

    def test_ensembl_gene(self):
        ''' Test ensembl GTF loading. '''
        call_command('pipeline', '--steps', 'stage', 'load', sections='ENSEMBL_GENE_GTF',
                     dir=TEST_DATA_DIR, ini=INI_FILE)
        idx = INI_CONFIG['ENSEMBL_GENE_GTF']['index']
        idx_type = INI_CONFIG['ENSEMBL_GENE_GTF']['index_type']
        elastic = Search(idx=idx, idx_type=idx_type)
        Search.index_refresh(idx)

        self.assertGreaterEqual(elastic.get_count()['count'], 1, "Count documents in the index")
        map_props = Gene.gene_mapping(idx, idx_type, test_mode=True).mapping_properties
        self._cmpMappings(elastic.get_mapping()[idx]['mappings'], map_props, idx_type)

    def test_gene_history_loader(self):
        ''' Test the gene history loading. '''
        call_command('pipeline', '--steps', 'load', sections='GENE_HISTORY',
                     dir=TEST_DATA_DIR, ini=INI_FILE)

        idx = INI_CONFIG['GENE_HISTORY']['index']
        idx_type = INI_CONFIG['GENE_HISTORY']['index_type']
        elastic = Search(idx=idx, idx_type=idx_type)
        Search.index_refresh(idx)

        self.assertTrue(elastic.get_count()['count'] > 1, "Count documents in the index")
        map_props = Gene.gene_history_mapping(idx, idx_type, test_mode=True).mapping_properties
        self._cmpMappings(elastic.get_mapping()[idx]['mappings'], map_props, idx_type)

    def _cmpMappings(self, map1, map2, idx_type):
        ''' Compare property keys and the types for two mappings. '''
        map1 = map1[idx_type]['properties']
        map2 = map2[idx_type]['properties']
        sorted_keys = list(map1.keys())
        sorted_keys.sort()
        for k in sorted_keys:
            self.assertTrue(k in map2, k)
            self.assertEqual(map1[k]['type'], map2[k]['type'], k)
