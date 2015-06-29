''' Tests for the download module. '''
from django.test import TestCase
from django.core.management import call_command
from data_pipeline.download import HTTPDownload, FTPDownload, MartDownload
from django.utils.six import StringIO
from elastic.elastic_settings import ElasticSettings
import os
import requests


class DownloadTest(TestCase):

    def test_ini_file(self):
        ''' Test ini file downloads. '''
        out = StringIO()
        ini_file = os.path.join(os.path.dirname(__file__), 'download.ini')
        call_command('pipeline', dir='/tmp', ini=ini_file, download=True, stdout=out)
        self.assertEqual(out.getvalue().strip(), "DOWNLOAD COMPLETE")

    def test_pub_ini_file(self):
        ''' Test publication ini file downloads. '''
        out = StringIO()
        ini_file = os.path.join(os.path.dirname(__file__), 'publication.ini')
        call_command('publications', '--dir', '/tmp', '--steps', 'download', 'load', ini=ini_file, stdout=out)
        self.assertEqual(out.getvalue().strip(), "DOWNLOAD COMPLETE")
        requests.delete(ElasticSettings.url() + '/' + 'test__publications')

    def test_file_cmd(self):
        out = StringIO()
        call_command('pipeline', dir='/tmp', url='http://t1dbase.org', download=True, stdout=out)
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
        call_command('pipeline', dir='/tmp', url='ftp://ftp.ebi.ac.uk/pub/databases/embl/README', download=True)
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
