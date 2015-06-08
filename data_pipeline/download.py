''' Download Data Module '''
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
import re

# Get an instance of a logger
logger = logging.getLogger(__name__)


def post_process(func):
    ''' Used as a decorator to apply L{PostProcess} functions. '''
    def wrapper(*args, **kwargs):
        success = func(*args, **kwargs)
        if success and 'section' in kwargs:
            section = kwargs['section']
            if 'post' in section:
                post_func = getattr(globals()['PostProcess'], section['post'])
                post_func(*args, **kwargs)
        return success
    return wrapper


class Download:
    ''' Handle data file downloads '''

    def download(self, url, dir_path, file_name=None, **kwargs):
        if file_name is None:
            file_name = self._url_to_file_name(url)

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

        config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        config.read(ini_file)

        success = False
        for section_name in config.sections():
            self._inherit_section(section_name, config)
            success = self._parse_ini(section_name, dir_path=dir_path, section=config[section_name])
        return success

    def _inherit_section(self, section_name, config):
        ''' Add in parameters from another config section when a double colon
        is found in the name. '''
        if '::' in section_name:
            inherit = section_name.split('::', maxsplit=1)[0]
            if inherit in config:
                for key in config[inherit]:
                    config[section_name][key] = config[inherit][key]

    @post_process
    def _parse_ini(self, fname, dir_path='.', section=None):
        success = False
        if 'output' in section:
            fname = section['output']

        if 'location' in section:
            if 'type' in section and section['type'] == 'emsembl_mart':
                qfilter = None
                if 'query_filter' in section:
                    qfilter = section['query_filter']
                elif 'ensgene_filter' in section:
                    qfilter = '<Filter name="ensembl_gene_id" value="%s"/>' % section['ensgene_filter']
                success = self.download(section['location'], dir_path, file_name=fname,
                                        tax=section['taxonomy'], attrs=section['attrs'],
                                        query_filter=qfilter, emsembl_mart=True)
            elif 'files' in section:
                files = section['files'].split(",")
                for f in files:
                    success = self.download(section['location']+"/"+f.strip(), dir_path)
            elif 'http_params' in section:
                success = self.download(section['location']+"?"+section['http_params'],
                                        dir_path, file_name=fname)
        return success

    def _url_to_file_name(self, url):
        name = url.split('/')[-1]
        if name == '':
            name = re.sub(r"[\/?\.:]", "", url)
        return name


class HTTPDownload:
    ''' HTTP downloader. '''

    @classmethod
    def download(cls, url, dir_path, file_name):
        r = requests.get(url, stream=True)
        if r.status_code != 200:
            logger.error("response "+str(r.status_code)+": "+url)
            return False

        monitor = Monitor(file_name, size=r.headers.get('content-length'))
        with open(os.path.join(dir_path, file_name), 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
                    monitor(chunk)
        return True


class FTPDownload:
    ''' FTP downloader. '''

    @classmethod
    def download(cls, url, dir_path, file_name):
        url_parse = urlparse(url)
        ftp_host = ftputil.FTPHost(url_parse.netloc, 'anonymous', '',
                                   session_factory=ftplib.FTP)
        size = ftp_host.path.getsize(url_parse.path)
        mon = Monitor(file_name, size=size)
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
            '<Query virtualSchemaName="default" formatter="TSV" ' \
            'header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">' \
            '<Dataset name="%s" interface="default">%s%s' \
            '</Dataset>' \
            '</Query>' % (url, tax, query_filter, attrs_str)
        return HTTPDownload.download(url_query, dir_path, file_name)


class PostProcess:

    @classmethod
    def zcat(cls, *args, **kwargs):
        ''' Combine a list of compressed files. '''
        section = kwargs['section']
        dir_path = kwargs['dir_path']
        out = os.path.join(kwargs['dir_path'], section['output'])

        files = section['files'].split(",")
        with open(out, 'wb') as outf:
            for fname in files:
                with open(os.path.join(dir_path, fname), 'rb') as infile:
                    for line in infile:
                        outf.write(line)
                os.remove(os.path.join(dir_path, fname))

    @classmethod
    def xmlparse(cls, *args, **kwargs):
        out = os.path.join(kwargs['dir_path'], kwargs['section']['output'])

        tree = ET.parse(out)
        idlist = tree.find("IdList")
        ids = list(idlist.iter("Id"))
        os.remove(out)
        with open(out, 'w') as outf:
            for i in ids:
                outf.write(i.text+'\n')


class Monitor:
    ''' Monitor download progress. '''

    def __init__(self, file_name, size=None):
        if size is not None:
            self.size = int(size)
        self.size_progress = 0
        self.previous = 0
        self.start = time.time()
        self.file_name = file_name
        print("%s" % file_name, end="", flush=True)

    def __call__(self, chunk):
        self.size_progress += len(chunk)

        if not hasattr(self, 'size'):
            print("\r[%s] %s" % (self.size_progress, self.file_name), end="", flush=True)
            return

        progress = int(self.size_progress/self.size * 100)
        if progress != self.previous and progress % 10 == 0:
            time_taken = time.time() - self.start
            eta = (time_taken / self.size_progress) * (self.size - self.size_progress)
            print("\r[%s%s] eta:%ss  %s  " % ('=' * int(progress/2),
                                              ' ' * (50-int(progress/2)),
                                              str(int(eta)), self.file_name), end="", flush=True)
            self.previous = progress
