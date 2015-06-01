from django.test import TestCase
from django.core.management import call_command
from data_pipeline.download import HTTPDownload, FTPDownload, MartDownload
import os


class DownloadTest(TestCase):

    ''' HTTP tests '''
    def test_file_cmd(self):
        call_command('download', dir='/tmp', url='http://t1dbase.org')

    def test_http(self):
        self.assertTrue(HTTPDownload.download('http://t1dbase.org', '/tmp', 't1d.tmp'),
                        'HTTP download test')

    ''' FTP tests '''
    def test_ftp_cmd(self):
        out = os.path.join('/tmp', 'README')
        if os.path.isfile(out):
            os.remove(out)
        call_command('download', dir='/tmp', url='ftp://ftp.ebi.ac.uk/pub/databases/embl/README')
        self.assertTrue(os.path.isfile(os.path.join('/tmp', 'README')), 'FTP test command')
        os.remove(out)

    def test_ftp(self):
        self.assertTrue(FTPDownload.download('ftp://ftp.ebi.ac.uk/pub/databases/embl/README',
                                             '/tmp', 'ftp.test'),
                        'FTP download test')

    ''' Mart tests '''
    def test_mart(self):
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