import requests
import ftputil
from urllib.parse import urlparse
import ftplib
import time
import csv
import configparser
import os
import xml.etree.ElementTree as ET
from builtins import classmethod
import logging
import sys

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

    def download(self, url, dir_path, file_name=None):
        if file_name is None:
            file_name = url.split('/')[-1]
        print(file_name)
        if url.startswith("ftp://"):
            success = FTPDownload.download(url, dir_path, file_name)
        else:
            success = HTTPDownload.download(url, dir_path, file_name)
        return success

    def download_ini(self, ini_file, dir_path):
        config = configparser.ConfigParser()
        config.read(ini_file)
        ini_path = os.path.dirname(os.path.abspath(ini_file))
        print(config.sections())
        for section_name in config.sections():
            self._process_section(section=config[section_name], name=section_name,
                                  dir_path=dir_path, ini_path=ini_path)

    @post_process
    def _process_section(self, section=None, name=None, dir_path='.', ini_path=None):
        if 'location' in section:
            if 'type' in section and section['type'] == 'emsembl_mart':
                MartDownload.download(section)
                sys.exit()
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

    def download_file(self, filename, dir_path):
        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) > 1:
                    if row[0].startswith("#"):
                        continue

                    if len(row) > 2:
                        files = row[2].split(",")
                        for f in files:
                            self.download(row[1]+"/"+f.strip(), dir_path)
                    else:
                        self.download(row[1], dir_path)


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

    @classmethod
    def download(cls, section):
        taxonomy = section['taxonomy']
        location = section['location']
        attrs = section['attrs']
        attrs_str = ''
        for attr in attrs.split(','):
            attrs_str += '''<Attribute name="%s"/>''' % attr.strip()

        gene = 'ENSG00000134242'
        query_filter = '<Filter name="ensembl_gene_id" value="%s"/>' % gene

        urlTemplate = \
            '%s?query=' \
            '<?xml version="1.0" encoding="UTF-8"?>' \
            '<!DOCTYPE Query>' \
            '<Query virtualSchemaName="default" formatter="CSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">' \
            '<Dataset name="%s" interface="default">%s%s' \
            '</Dataset>' \
            '</Query>'
        martURL = urlTemplate % (location, taxonomy, query_filter, attrs_str)
        r = requests.get(martURL, stream=True)
        for line in r.iter_lines():
            print (line)


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

    def __init__(self, size):
        if size is not None:
            self.size = int(size)
        else:
            self.size = 1000000
        self.size_progress = 0
        self.previous = 0
        self.start = time.time()

    def __call__(self, chunk):
        self.size_progress += len(chunk)
        percent_progress = int(self.size_progress/self.size * 100)
        if percent_progress != self.previous and percent_progress % 10 == 0:
            time_taken = time.time() - self.start
            eta = (time_taken / self.size_progress) * (self.size - self.size_progress)
            print("\r[%s%s] eta:%ss    " % ('=' * int(percent_progress/2),
                                            ' ' * (50-int(percent_progress/2)),
                                            str(int(eta))), end="")
            self.previous = percent_progress
