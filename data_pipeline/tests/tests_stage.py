''' Tests for the staging. '''
from django.test import TestCase
from django.utils.six import StringIO
import os
import data_pipeline
import shutil
import logging
logger = logging.getLogger(__name__)


def setUpModule():
    logging.debug('Module setup...')
    pass


def tearDownModule():
    logging.debug('Module teardown...')
    app_data_dir = os.path.dirname(data_pipeline.__file__)
    test_data_dir = app_data_dir + '/tests/data'
    stage_data_dir = test_data_dir + '/STAGE'
    if(os.path.exists(stage_data_dir)):
        shutil.rmtree(stage_data_dir)


class StagingTest(TestCase):

    def setUp(self):
        logging.debug('Class setup...' + self.__class__.__name__)
        print('Class setup...' + self.__class__.__name__)
        self.out = StringIO()
        self.ini_file = os.path.join(os.path.dirname(__file__), 'download.ini')
        self.app_data_dir = os.path.dirname(data_pipeline.__file__)
        self.test_data_dir = self.app_data_dir + '/tests/data'
        logging.debug(self.test_data_dir)

    '''Add some generic staging tests here'''
    def test_dummy(self):
        self.assertTrue(1 < 2, "Dummy Test")
