''' Tests for the staging and loading. '''
from django.test import TestCase
from django.core.management import call_command
import os
import data_pipeline
import shutil
from data_pipeline.utils import IniParser
from elastic.elastic_settings import ElasticSettings
import requests
from elastic.search import Search, ElasticQuery
from data_pipeline.helper.gene import Gene
import logging
import json
from elastic.query import Query, TermsFilter

# Get an instance of a logger
logger = logging.getLogger(__name__)
IDX_SUFFIX = ElasticSettings.getattr('TEST')
MY_INI_FILE = os.path.join(os.path.dirname(__file__), IDX_SUFFIX + '_test_download.ini')
TEST_DATA_DIR = os.path.dirname(data_pipeline.__file__) + '/tests/data'


def setUpModule():
    ''' Change ini config (MY_INI_FILE) to use the test suffix when
    creating pipeline indices. '''
    ini_file = os.path.join(os.path.dirname(__file__), 'test_download.ini')
    if os.path.isfile(MY_INI_FILE):
        return

    with open(MY_INI_FILE, 'w') as new_file:
        with open(ini_file) as old_file:
            for line in old_file:
                new_file.write(line.replace('auto_tests', IDX_SUFFIX))


def tearDownModule():
    if os.path.exists(TEST_DATA_DIR + '/STAGE'):
        shutil.rmtree(TEST_DATA_DIR + '/STAGE')
    # remove index created
    INI_CONFIG = IniParser().read_ini(MY_INI_FILE)
    requests.delete(ElasticSettings.url() + '/' + INI_CONFIG['GENE_HISTORY']['index'])
    requests.delete(ElasticSettings.url() + '/' + INI_CONFIG['DBSNP']['index'])
    os.remove(MY_INI_FILE)
    ens_dir = os.path.join(TEST_DATA_DIR, 'DOWNLOAD', 'ENSMART_GENE')
    if os.path.exists(ens_dir):
        shutil.rmtree(ens_dir)


class StageTest(TestCase):

    def test_gene_interaction_parse(self):
        ''' Test ini file downloads. '''
        call_command('pipeline', '--steps', 'stage', sections='INTACT',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)

        self.assertTrue(os.path.exists(TEST_DATA_DIR + '/STAGE/INTACT'))
        self.assertTrue(os.path.isfile(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.out'))
        self.assertTrue(os.path.isfile(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.json'))
        self.assertTrue(os.access(TEST_DATA_DIR + '/STAGE/INTACT/intact.zip.json', os.R_OK))


class LoadTest(TestCase):

    def test_gene_pipeline(self):
        ''' Test gene pipeline. '''

        INI_CONFIG = IniParser().read_ini(MY_INI_FILE)
        idx = INI_CONFIG['ENSEMBL_GENE_GTF']['index']
        idx_type = INI_CONFIG['ENSEMBL_GENE_GTF']['index_type']

        ''' 1. Test ensembl GTF loading. '''
        call_command('pipeline', '--steps', 'stage', 'load', sections='ENSEMBL_GENE_GTF',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)

        elastic = Search(idx=idx, idx_type=idx_type)
        self.assertGreaterEqual(elastic.get_count()['count'], 1, "Count documents in the index")
        map1_props = Gene.gene_mapping(idx, idx_type, test_mode=True).mapping_properties
        map2_props = elastic.get_mapping()
        if idx not in map2_props:
            logger.error("MAPPING ERROR: "+json.dumps(map2_props))
        self._cmpMappings(map2_props[idx]['mappings'], map1_props, idx_type)

        ''' 2. Test adding entrez ID to documents '''
        call_command('pipeline', '--steps', 'load', sections='GENE2ENSEMBL',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)
        query = ElasticQuery.query_string("PTPN22", fields=["symbol"])
        elastic = Search(query, idx=idx)
        docs = elastic.search().docs
        self.assertEqual(len(docs), 1)
        self.assertTrue('entrez' in getattr(docs[0], "dbxrefs"))
        self.assertEqual(getattr(docs[0], "dbxrefs")['entrez'], '26191')

        ''' 3. Add uniprot and fill in missing entrez fields. '''
        call_command('pipeline', '--steps', 'download', 'load', sections='ENSMART_GENE',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)
        query = ElasticQuery.query_string("DNMT3L", fields=["symbol"])
        elastic = Search(query, idx=idx)
        docs = elastic.search().docs
        self.assertTrue('entrez' in getattr(docs[0], "dbxrefs"))
        self.assertTrue('swissprot' in getattr(docs[0], "dbxrefs"))

        ''' 4. Add gene synonyms and dbxrefs. '''
        call_command('pipeline', '--steps', 'load', sections='GENE_INFO',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)
        query = ElasticQuery.query_string("PTPN22", fields=["symbol"])
        elastic = Search(query, idx=idx)
        docs = elastic.search().docs
        self.assertTrue('PTPN8' in getattr(docs[0], "synonyms"))

        ''' 5. Add PMIDs to gene docs. '''
        call_command('pipeline', '--steps', 'load', sections='GENE_PUBS',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)
        query = ElasticQuery.query_string("PTPN22", fields=["symbol"])
        elastic = Search(query, idx=idx)
        docs = elastic.search().docs
        self.assertGreater(len(getattr(docs[0], "pmids")), 0)

        ''' 6. Add ortholog data. '''
        call_command('pipeline', '--steps', 'load', sections='ENSMART_HOMOLOG',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)
        query = ElasticQuery.query_string("PTPN22", fields=["symbol"])
        elastic = Search(query, idx=idx)
        docs = elastic.search().docs
        dbxrefs = getattr(docs[0], "dbxrefs")
        self.assertTrue('orthologs' in dbxrefs, dbxrefs)
        self.assertTrue('mmusculus' in dbxrefs['orthologs'], dbxrefs)
        self.assertEqual('ENSMUSG00000027843', dbxrefs['orthologs']['mmusculus']['ensembl'])

        query = ElasticQuery.filtered(Query.match_all(),
                                      TermsFilter.get_terms_filter("dbxrefs.orthologs.mmusculus.ensembl",
                                                                   ['ENSMUSG00000027843']))
        docs = Search(query, idx=idx, size=1).search().docs
        self.assertEqual(len(docs), 1)

        ''' 7. Add mouse ortholog link to MGI '''
        call_command('pipeline', '--steps', 'load', sections='ENSEMBL2MGI',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        Search.index_refresh(idx)
        docs = Search(query, idx=idx, size=1).search().docs
        dbxrefs = getattr(docs[0], "dbxrefs")
        self.assertEqual('ENSMUSG00000027843', dbxrefs['orthologs']['mmusculus']['ensembl'])
        self.assertEqual('107170', dbxrefs['orthologs']['mmusculus']['MGI'])

    def test_marker_pipeline(self):
        ''' Test marker pipeline. '''
        call_command('pipeline', '--steps', 'load', sections='DBSNP',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)

        INI_CONFIG = IniParser().read_ini(MY_INI_FILE)
        idx = INI_CONFIG['DBSNP']['index']
        idx_type = INI_CONFIG['DBSNP']['index_type']
        elastic = Search(idx=idx, idx_type=idx_type)
        Search.index_refresh(idx)
        self.assertGreater(elastic.get_count()['count'], 0)

        call_command('pipeline', '--steps', 'load', sections='RSMERGEARCH',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)
        idx = INI_CONFIG['RSMERGEARCH']['index']
        idx_type = INI_CONFIG['RSMERGEARCH']['index_type']
        elastic = Search(idx=idx, idx_type=idx_type)
        Search.index_refresh(idx)
        self.assertGreater(elastic.get_count()['count'], 0)

    def test_gene_history_loader(self):
        ''' Test the gene history loading. '''
        call_command('pipeline', '--steps', 'load', sections='GENE_HISTORY',
                     dir=TEST_DATA_DIR, ini=MY_INI_FILE)

        INI_CONFIG = IniParser().read_ini(MY_INI_FILE)
        idx = INI_CONFIG['GENE_HISTORY']['index']
        idx_type = INI_CONFIG['GENE_HISTORY']['index_type']
        elastic = Search(idx=idx, idx_type=idx_type)
        Search.index_refresh(idx)

        self.assertTrue(elastic.get_count()['count'] > 1, "Count documents in the index")
        map1_props = Gene.gene_history_mapping(idx, idx_type, test_mode=True).mapping_properties
        map2_props = elastic.get_mapping()
        if idx not in map2_props:
            logger.error("MAPPING ERROR: "+json.dumps(map2_props))
        self._cmpMappings(map2_props[idx]['mappings'], map1_props, idx_type)

    def _cmpMappings(self, map1, map2, idx_type):
        ''' Compare property keys and the types for two mappings. '''
        map1 = map1[idx_type]['properties']
        map2 = map2[idx_type]['properties']
        sorted_keys = list(map1.keys())
        sorted_keys.sort()
        for k in sorted_keys:
            self.assertTrue(k in map2, k)
            if 'type' in map1[k]:
                self.assertEqual(map1[k]['type'], map2[k]['type'], k)
