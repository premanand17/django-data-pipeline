''' Tests for the staging. '''
from django.test import TestCase
from django.core.management import call_command
from django.utils.six import StringIO
import os
import data_pipeline
import shutil
import logging
logger = logging.getLogger(__name__)


def tearDownModule():
    logging.debug('Module teardown...')
    app_data_dir = os.path.dirname(data_pipeline.__file__)
    test_data_dir = app_data_dir + '/tests/data'
    shutil.rmtree(test_data_dir + '/STAGE')


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
    pass


class GenePathwayProcessTest(TestCase):
    pass




