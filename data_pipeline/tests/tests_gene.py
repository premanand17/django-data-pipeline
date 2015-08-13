''' Tests for the staging. '''
from django.test import TestCase
from django.core.management import call_command
from django.utils.six import StringIO
import os
import data_pipeline
import shutil
import logging
from data_pipeline.helper.gene_interactions import GeneInteractions
import pprint
import json
logger = logging.getLogger(__name__)


def tearDownModule():
    logging.debug('Module teardown...')
    app_data_dir = os.path.dirname(data_pipeline.__file__)
    test_data_dir = app_data_dir + '/tests/data'
    stage_data_dir = test_data_dir + '/STAGE'
    #if(os.path.exists(stage_data_dir)):
    #    shutil.rmtree(stage_data_dir)


class GeneInteractionStagingTest(TestCase):

    def setUp(self):
        logging.debug('Class setup...' + self.__class__.__name__)
        self.out = StringIO()
        self.ini_file = os.path.join(os.path.dirname(__file__), 'download.ini')
        self.app_data_dir = os.path.dirname(data_pipeline.__file__)
        self.test_data_dir = self.app_data_dir + '/tests/data'
        logging.debug(self.test_data_dir)

    def test_intact(self):
        ''' Test intact staging. '''
        logging.debug("Test intact staging")
        call_command('pipeline', '--steps', 'stage', sections='INTACT',
                     dir=self.test_data_dir, ini=self.ini_file, stdout=self.out)

        self.assertTrue(os.path.exists(self.test_data_dir + '/STAGE/INTACT'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/INTACT/intact.zip.out'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/INTACT/intact.zip.json'))
        self.assertTrue(os.access(self.test_data_dir + '/STAGE/INTACT/intact.zip.json', os.R_OK))

        # check if the .out file is in right format - this file will have interactors with ensemblid
        with open(self.test_data_dir + '/STAGE/INTACT/intact.zip.out', 'r') as f:
            header = f.readline()
            first_line = f.readline()
            f.close()

        self.assertEqual(header, "interactorA\tinteractorB\tpubmed\n", "header is ok")
        self.assertEqual(first_line, "ENSG00000078053\tENSG00000159082\t10542231\n", "first line is ok")

        # check if the json outputput file has the right data
        with open(self.test_data_dir + '/STAGE/INTACT/intact.zip.json', 'r') as f:
            docs_line = f.readline()
            interactor_line = f.readline()
            f.close()

        self.assertIn("docs", docs_line, "docs line OK")
        self.assertIn("interaction_source", interactor_line, "contains interaction_source")
        self.assertIn("intact", interactor_line, "contains intact as interaction_source")

    def test_bioplex(self):
        ''' Test bioplex staging. '''
        logging.debug("Test bioplex staging")
        call_command('pipeline', '--steps', 'stage', sections='BIOPLEX',
                     dir=self.test_data_dir, ini=self.ini_file, stdout=self.out)

        self.assertTrue(os.path.exists(self.test_data_dir + '/STAGE/BIOPLEX'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.out'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.json'))
        self.assertTrue(os.access(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.json', os.R_OK))

        # check if the .out file is in right format - this file will have interactors with ensemblid
        with open(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.out', 'r') as f:
            header = f.readline()
            first_line = f.readline()
            f.close()

        self.assertEqual(header, "interactorA\tinteractorB\n", "header is ok")
        self.assertEqual(first_line, "ENSG00000196839\tENSG00000196604\n", "first line is ok")

        # check if the json outputput file has the right data
        with open(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.json', 'r') as f:
            docs_line = f.readline()
            interactor_line = f.readline()
            f.close()

        self.assertIn("docs", docs_line, "docs line OK")
        self.assertIn("interaction_source", interactor_line, "contains interaction_source")
        self.assertIn("bioplex", interactor_line, "contains bioplex as interaction_source")


class GenePathwayStagingTest(TestCase):

    def setUp(self):
        logging.debug('Class setup...' + self.__class__.__name__)
        self.out = StringIO()
        self.ini_file = os.path.join(os.path.dirname(__file__), 'download.ini')
        self.app_data_dir = os.path.dirname(data_pipeline.__file__)
        self.test_data_dir = self.app_data_dir + '/tests/data'
        logging.debug(self.test_data_dir)

    def test_msigdb(self):
        ''' Test msigdb staging. '''
        logging.debug("Test msigdb staging")
        call_command('pipeline', '--steps', 'stage', sections='MSIGDB',
                     dir=self.test_data_dir, ini=self.ini_file, stdout=self.out)

        self.assertTrue(os.path.exists(self.test_data_dir + '/STAGE/MSIGDB'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.biocarta.v5.0.entrez.gmt.json'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.kegg.v5.0.entrez.gmt.json'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.reactome.v5.0.entrez.gmt.json'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c5.all.v5.0.entrez.gmt.json'))

        self.assertTrue(os.access(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.biocarta.v5.0.entrez.gmt.json', os.R_OK))


class GeneInteractionProcessTest(TestCase):
    '''Test functions in GeneInteractions class'''

    def test_group_binary_interactions(self):
        '''
        Test if the interactors are grouped correctly
        '''
        binary_interactions = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (3, 4),
                               (4, 5), (5, 6), (5, 7), (6, 8), (7, 8), (8, 9)]

        expected_dict = {
            '0': ['1', '2'],
            '1': ['0', '2', '3'],
            '2': ['0', '1', '3'],
            '3': ['1', '2', '4'],
            '4': ['3', '5'],
            '5': ['4', '6', '7'],
            '6': ['5', '8'],
            '7': ['5', '8'],
            '8': ['6', '7', '9'],
            '9': ['8']
            }

        grouped_interactions = GeneInteractions._group_binary_interactions(binary_interactions)
        self.assertDictEqual(grouped_interactions, expected_dict, "Grouping of interactions done OK")
        self.assertTrue(10 == len(grouped_interactions), "dict contains 10 groups")
        self.assertIn('8', expected_dict, "key 8 is in dict")
        self.assertIn('9', grouped_interactions['8'], "value 9 is in dict['8'] ")
        self.assertIn('1', grouped_interactions['0'], "value 1 is in dict['0'] ")

    def test_interactor_json_decorator(self):
        json_a = GeneInteractions.interactor_json_decorator("geneA")
        self.assertDictEqual(json_a, {'interactor': 'geneA'}, "JSON equal for geneA")

    def test_interaction_json_decorator(self):
        json_interaction = GeneInteractions.interaction_json_decorator("intact", "geneA", ["geneB", "geneC", "geneD"],
                                                                       "pubmed", ["1234", "4567", "8912"])

        self.assertIn('"_parent": "geneA"', json_interaction, "_parent found in json")
        self.assertIn('"interaction_source": "intact"', json_interaction, "_parent found in json")

        json_interaction = GeneInteractions.interaction_json_decorator("bioplex", "geneA", ["geneB", "geneC", "geneD"])
        self.assertJSONEqual(json_interaction,
                             json.dumps({"interactors": [{"interactor": "geneB"},
                                                         {"interactor": "geneC"},
                                                         {"interactor": "geneD"}],
                                         "interaction_source": "bioplex", "_parent": "geneA"}),
                             "JSON for interaction equal")


class GenePathwayProcessTest(TestCase):
    pass # ToDO
