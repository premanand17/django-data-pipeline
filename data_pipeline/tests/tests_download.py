''' Tests for the download module. '''
from django.test import TestCase
from django.core.management import call_command
from data_pipeline.download import HTTPDownload, FTPDownload, MartDownload
from django.utils.six import StringIO
from elastic.elastic_settings import ElasticSettings
import os
import requests
from data_pipeline.utils import IniParser

IDX_SUFFIX = ElasticSettings.getattr('TEST')
MY_PUB_INI_FILE = os.path.join(os.path.dirname(__file__), IDX_SUFFIX + '_test_publication.ini')


def setUpModule():
    ''' Change ini config (MY_PUB_INI_FILE) to use the test suffix when
    creating publication pipeline index. '''
    ini_file = os.path.join(os.path.dirname(__file__), 'test_publication.ini')
    if os.path.isfile(MY_PUB_INI_FILE):
        return

    with open(MY_PUB_INI_FILE, 'w') as new_file:
        with open(ini_file) as old_file:
            for line in old_file:
                new_file.write(line.replace('auto_tests', IDX_SUFFIX))


def tearDownModule():
    # remove index created
    INI_CONFIG = IniParser().read_ini(MY_PUB_INI_FILE)
    requests.delete(ElasticSettings.url() + '/' + INI_CONFIG['DISEASE']['index'])
    os.remove(MY_PUB_INI_FILE)


class DownloadTest(TestCase):

    def test_ini_file(self):
        ''' Test ini file downloads. '''
        out = StringIO()
        ini_file = os.path.join(os.path.dirname(__file__), 'test_download.ini')
        call_command('pipeline', '--steps', 'download', sections='ENSEMBL_GENE', dir='/tmp', ini=ini_file, stdout=out)
        self.assertEqual(out.getvalue().strip(), "DOWNLOAD COMPLETE")

    def test_pub_ini_file(self):
        ''' Test publication ini file downloads. '''
        out = StringIO()
        call_command('publications', '--dir', '/tmp', '--steps', 'download', 'load', ini=MY_PUB_INI_FILE, stdout=out)
        self.assertEqual(out.getvalue().strip(), "DOWNLOAD COMPLETE")

    def test_file_cmd(self):
        out = StringIO()
        call_command('pipeline', '--steps', 'download', dir='/tmp', url='http://t1dbase.org', stdout=out)
        self.assertEqual(out.getvalue().strip(), "DOWNLOAD COMPLETE")

    def test_http(self):
        ''' Test downloading over HTTP. '''
        self.assertTrue(HTTPDownload.download('http://t1dbase.org', '/tmp', 't1d.tmp'),
                        'HTTP download test')

    def test_ftp_cmd(self):
        ''' Test downloading over FTP. '''
        out = os.path.join('/tmp', 'README')
        if os.path.isfile(out):
            os.remove(out)
        call_command('pipeline', '--steps', 'download', dir='/tmp', url='ftp://ftp.ebi.ac.uk/pub/databases/embl/README')
        self.assertTrue(os.path.isfile(os.path.join('/tmp', 'README')), 'FTP test command')
        os.remove(out)

    def test_ftp_exists(self):
        ''' Test FTP exists. '''
        self.assertTrue(FTPDownload.exists('ftp://ftp.ebi.ac.uk/'), 'FTP file/dir exists')
        self.assertFalse(FTPDownload.exists('ftp://ftp.ebi.ac.uk/xxxx'), 'FTP file/dir exists')

    def test_ftp_mtime(self):
        ''' Test mtime from a file on a FTP server. '''
        self.assertTrue(FTPDownload.mtime('ftp://ftp.ebi.ac.uk/pub/databases/embl/README') > 0,
                        'FTP file/dir exists')

    def test_ftp(self):
        ''' Test downloading over FTP. '''
        self.assertTrue(FTPDownload.download('ftp://ftp.ebi.ac.uk/pub/databases/embl/README',
                                             '/tmp', 'ftp.test'),
                        'FTP download test')

    def test_mart(self):
        ''' Test downloading from MART. '''
        query_filter = '<Filter name="ensembl_gene_id" value="ENSG00000134242"/>'
        attrs = '<Attribute name="ensembl_gene_id"/>' \
                '<Attribute name="hgnc_symbol"/>' \
                '<Attribute name="gene_biotype"/>' \
                '<Attribute name="start_position"/>' \
                '<Attribute name="end_position"/>' \
                '<Attribute name="strand"/>'

        self.assertTrue(
            MartDownload.download('http://ensembl.org/biomart/martservice/', '/tmp',
                                  'hsapiens_gene_ensembl.out',
                                  query_filter=query_filter,
                                  tax='hsapiens_gene_ensembl', attrs=attrs),
            'Mart download')
