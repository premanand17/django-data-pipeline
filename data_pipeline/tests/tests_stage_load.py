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


class StagingTest(TestCase):

    def test_gene_interaction_parse(self):
        ''' Test ini file downloads. '''
        call_command('pipeline', '--steps', 'stage', sections='INTACT',
                     dir=TEST_DATA_DIR, ini=INI_FILE)

        self.assertTrue(os.path.exists(TEST_DATA_DIR + '/STAGE/INTACT'))
        self.assertTrue(os.path.isfile(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.out'))
        self.assertTrue(os.path.isfile(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.json'))
        self.assertTrue(os.access(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.json', os.R_OK))

    def test_gene_history_loader(self):
        call_command('pipeline', '--steps', 'load', sections='GENE_HISTORY',
                     dir=TEST_DATA_DIR, ini=INI_FILE)

        idx = INI_CONFIG['GENE_HISTORY']['index']
        idx_type = INI_CONFIG['GENE_HISTORY']['index_type']
        elastic = Search(idx=idx)
        Search.index_refresh(idx)
        mapping = elastic.get_mapping()

        self.assertTrue(elastic.get_count()['count'] > 1, "Count documents in the index")
        map_props = Gene.gene_history_mapping(idx, idx_type, test_mode=True).mapping_properties[idx_type]['properties']
        load_map_props = mapping[idx]['mappings'][idx_type]['properties']

        sorted_keys = list(load_map_props.keys())
        sorted_keys.sort()
        for k in sorted_keys:
            self.assertTrue(k in map_props)
            self.assertEqual(load_map_props[k]['type'], map_props[k]['type'])
