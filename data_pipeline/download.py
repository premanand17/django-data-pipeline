''' Download Data Module '''
from urllib.parse import urlparse
from builtins import classmethod
import requests
import ftputil
import ftplib
import os
import logging
import re
from .utils import IniParser
from .utils import post_process
from .utils import Monitor

# Get an instance of a logger
logger = logging.getLogger(__name__)


class Download(IniParser):
    ''' Handle data file downloads '''

    def download(self, url, dir_path, file_name=None, **kwargs):
        if file_name is None:
            file_name = self._url_to_file_name(url)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        if url.startswith("ftp://"):
            success = FTPDownload.download(url, dir_path, file_name, **kwargs)
        elif 'emsembl_mart' in kwargs:
            success = MartDownload.download(url, dir_path, file_name, **kwargs)
        else:
            success = HTTPDownload.download(url, dir_path, file_name, **kwargs)
        print()
        return success

    def download_ini(self, ini_file, dir_path, sections=None):
        ''' Download data defined in the ini file. '''
        return self.process_sections(self.read_ini(ini_file), dir_path, sections)

    @post_process
    def process_section(self, fname, section_dir_name, base_dir_path,
                        dir_path='.', section=None, stage='download'):
        ''' Overrides L{IniParser.process_section} to process a section
        in the config file '''
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


class HTTPDownload(object):
    ''' HTTP downloader. '''

    @classmethod
    def download(cls, url, dir_path, file_name, append=False):
        r = requests.get(url, stream=True, timeout=10)
        if r.status_code != 200:
            logger.error("response "+str(r.status_code)+": "+url)
            return False

        monitor = Monitor(file_name, size=r.headers.get('content-length'))
        if append:
            access = 'ab'
        else:
            access = 'wb'

        with open(os.path.join(dir_path, file_name), access) as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()
                    monitor(chunk)
        r.close()
        return True

    @classmethod
    def status(cls, url):
        return requests.get(url).status_code


class FTPDownload(object):
    ''' FTP downloader. '''

    @classmethod
    def download(cls, url, dir_path, file_name, username='anonymous', password=''):
        url_parse = urlparse(url)
        ftp_host = ftputil.FTPHost(url_parse.netloc, username, password,
                                   session_factory=ftplib.FTP)
        size = ftp_host.path.getsize(url_parse.path)
        mon = Monitor(file_name, size=size)
        ftp_host.download(url_parse.path, os.path.join(dir_path, file_name), callback=mon)
        ftp_host.close()

        if mon.size_progress != size:
            logger.error(file_name)
            logger.error("download size: "+mon.size_progress+" server size: "+size)
        return mon.size_progress == size

    @classmethod
    def mtime(cls, url, username='anonymous', password=''):
        ''' Time of most recent content modification in seconds '''
        url_parse = urlparse(url)
        ftp_host = ftputil.FTPHost(url_parse.netloc, username, password,
                                   session_factory=ftplib.FTP)
        return getattr(ftp_host.stat(url_parse.path), 'st_mtime')

    @classmethod
    def exists(cls, url, username='anonymous', password=''):
        url_parse = urlparse(url)
        ftp_host = ftputil.FTPHost(url_parse.netloc, username, password,
                                   session_factory=ftplib.FTP)
        return ftp_host.path.exists(url_parse.path)


class MartDownload(object):
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
        @type  qobjectuery_filter: string
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
            'header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">' \
            '<Dataset name="%s" interface="default">%s%s' \
            '</Dataset>' \
            '</Query>' % (url, tax, query_filter, attrs_str)
        return HTTPDownload.download(url_query, dir_path, file_name)
