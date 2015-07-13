''' Tests for the staging. '''
from django.test import TestCase
from django.core.management import call_command
from django.utils.six import StringIO
import os
import data_pipeline
import shutil


class StagingTest(TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        app_data_dir = os.path.dirname(data_pipeline.__file__)
        test_data_dir = app_data_dir + '/tests/data'
        shutil.rmtree(test_data_dir + '/STAGE')

    def test_gene_interaction_parse(self):
        ''' Test ini file downloads. '''
        out = StringIO()
        ini_file = os.path.join(os.path.dirname(__file__), 'download.ini')
        app_data_dir = os.path.dirname(data_pipeline.__file__)
        test_data_dir = app_data_dir + '/tests/data'
        print(test_data_dir)

        call_command('pipeline', '--steps', 'stage', sections='INTACT', dir=test_data_dir, ini=ini_file, stdout=out)

        self.assertTrue(os.path.exists(test_data_dir + '/STAGE/INTACT'))
        self.assertTrue(os.path.isfile(test_data_dir + '/STAGE/INTACT/intact.zip.out'))
        self.assertTrue(os.path.isfile(test_data_dir + '/STAGE/INTACT/intact.zip.json'))
        self.assertTrue(os.access(test_data_dir + '/STAGE/INTACT/intact.zip.json', os.R_OK))
