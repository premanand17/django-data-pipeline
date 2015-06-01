from urllib.parse import urlparse
from builtins import classmethod
import xml.etree.ElementTree as ET
import requests
import ftputil
import ftplib
import time
import configparser
import os
import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)


def post_process(func):
    def inner(*args, **kwargs):
        success = func(*args, **kwargs)
        if success and 'section' in kwargs:
            section = kwargs['section']
            if 'post' in section:
                post_func = getattr(globals()['PostProcess'], section['post'])
                post_func(*args, **kwargs)
    return inner


class Download:
    ''' Handle data file downloads '''

    def download(self, url, dir_path, file_name=None, **kwargs):
        if file_name is None:
            file_name = url.split('/')[-1]
        print(file_name)
        if url.startswith("ftp://"):
            success = FTPDownload.download(url, dir_path, file_name)
        elif 'emsembl_mart' in kwargs:
            success = MartDownload.download(url, dir_path, file_name, **kwargs)
        else:
            success = HTTPDownload.download(url, dir_path, file_name)
        print()
        return success

    def download_ini(self, ini_file, dir_path):
        ''' Download data defined in the ini file. '''

        # check for ini file in data pipeline home
        if not os.path.isfile(ini_file):
            DOWNLOAD_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
            tmp = os.path.join(DOWNLOAD_BASE_DIR, 'data_pipeline', ini_file)
            if os.path.isfile(tmp):
                ini_file = tmp

        config = configparser.ConfigParser()
        config.read(ini_file)
        ini_path = os.path.dirname(os.path.abspath(ini_file))
        print(config.sections())
        for section_name in config.sections():
            self._process_ini_section(section=config[section_name], name=section_name,
                                      dir_path=dir_path, ini_path=ini_path)

    @post_process
    def _process_ini_section(self, section=None, name=None, dir_path='.', ini_path=None):
        if 'location' in section:
            if 'type' in section and section['type'] == 'emsembl_mart':
                file_name = name
                if 'output' in section:
                    file_name = section['output']
                query_filter = None
                if 'query_filter' in section:
                    query_filter = section['query_filter']
                elif 'ensgene_filter' in section:
                    query_filter = '<Filter name="ensembl_gene_id" value="%s"/>' % section['ensgene_filter']
                self.download(section['location'], dir_path, file_name=file_name,
                              tax=section['taxonomy'], attrs=section['attrs'],
                              query_filter=query_filter, emsembl_mart=True)
            elif 'files' in section:
                files = section['files'].split(",")
                for f in files:
                    self.download(section['location']+"/"+f.strip(), dir_path)
            elif 'http_params' in section:
                file_name = name
                if 'output' in section:
                    file_name = section['output']
                self.download(section['location']+"?"+section['http_params'],
                              dir_path, file_name=file_name)


class HTTPDownload:

    @classmethod
    def download(cls, url, dir_path, file_name):
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            logger.error("response "+str(r.status_code)+": "+url)
            return False

        monitor = Monitor(r.headers.get('content-length'))
        with open(os.path.join(dir_path, file_name), 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
                    monitor(chunk)
        return True


class FTPDownload:

    @classmethod
    def download(cls, url, dir_path, file_name):
        url_parse = urlparse(url)
        ftp_host = ftputil.FTPHost(url_parse.netloc, 'anonymous', '',
                                   session_factory=ftplib.FTP)
        size = ftp_host.path.getsize(url_parse.path)
        mon = Monitor(size)
        ftp_host.download(url_parse.path, os.path.join(dir_path, file_name), callback=mon)

        if mon.size_progress != size:
            logger.error(file_name)
            logger.error("download size: "+mon.size_progress+" server size: "+size)
        return mon.size_progress == size


class MartDownload:
    ''' Biomart webservice downloads. '''

    @classmethod
    def download(cls, url, dir_path, file_name,
                 query_filter='', tax='', attrs='', **kwargs):
        '''
        @type  url: str
        @param url: The location of the mart service.
        @type  dir_path: str
        @param dir_path: Directory to write results to.
        @type  file_name: integer
        @param file_name: Output file name.
        @type  query_filter: string
        @keyword query_filter: Filter to be applied
                  (e.g. <Filter name="ensembl_gene_id" value="ENSG00000134242"/>.
        @type  tax: string
        @keyword tax: Taxonomy
        @type  attrs: string
        @keyword attrs: Comma separated attributes
        '''
        attrs_str = ''.join('<Attribute name="%s"/>' % a.strip() for a in attrs.split(','))
        url_query = \
            '%s?query=' \
            '<?xml version="1.0" encoding="UTF-8"?>' \
            '<!DOCTYPE Query>' \
            '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">' \
            '<Dataset name="%s" interface="default">%s%s' \
            '</Dataset>' \
            '</Query>' % (url, tax, query_filter, attrs_str)
        return HTTPDownload.download(url_query, dir_path, file_name)


class PostProcess:

    @classmethod
    def zcat(cls, *args, **kwargs):
        ''' Combine a list of compressed files. '''
        section = kwargs['section']
        outfile = section['output']
        dir_path = kwargs['dir_path']
        files = section['files'].split(",")
        with open(dir_path+"/"+outfile, 'wb') as outf:
            for fname in files:
                with open(dir_path+"/"+fname, 'rb') as infile:
                    for line in infile:
                        outf.write(line)
            os.remove(dir_path+"/"+fname)

    @classmethod
    def xmlparse(cls, *args, **kwargs):
        section = kwargs['section']
        outfile = section['output']
        dir_path = kwargs['dir_path']

        tree = ET.parse(dir_path+"/"+outfile)
        idlist = tree.find("IdList")
        ids = list(idlist.iter("Id"))
        os.remove(dir_path+"/"+outfile)
        with open(dir_path+"/"+outfile, 'w') as outf:
            for i in ids:
                outf.write(i.text+'\n')


class Monitor:
    ''' Monitor download progress. '''

    def __init__(self, size=None):
        if size is not None:
            self.size = int(size)
        self.size_progress = 0
        self.previous = 0
        self.start = time.time()

    def __call__(self, chunk):
        self.size_progress += len(chunk)

        if not hasattr(self, 'size'):
            print("\r[%s]" % self.size_progress, end="")
            return

        percent_progress = int(self.size_progress/self.size * 100)
        if percent_progress != self.previous and percent_progress % 10 == 0:
            time_taken = time.time() - self.start
            eta = (time_taken / self.size_progress) * (self.size - self.size_progress)
            print("\r[%s%s] eta:%ss    " % ('=' * int(percent_progress/2),
                                            ' ' * (50-int(percent_progress/2)),
                                            str(int(eta))), end="")
            self.previous = percent_progress