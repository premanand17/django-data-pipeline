''' Helper functions and classes '''
import os
import configparser
import time
import xml.etree.ElementTree as ET

from elastic.search import Search, ElasticQuery
from elastic.query import Query
from .helper.pubs import Pubs
import json
from elastic.management.loaders.loader import Loader


class Monitor(object):
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


def process_wrapper(*args, **kwargs):
    section = kwargs['section']
    ini_tag = None
    if kwargs['stage'] == 'Download':
        ini_tag = 'post'
    elif kwargs['stage'] == 'Stage':
        ini_tag = 'stage'
    elif 'Load' in kwargs['stage']:
        ini_tag = 'load'
    if ini_tag is not None:
        if ini_tag in section:
            post_func = getattr(globals()['PostProcess'], section[ini_tag])
            post_func(*args, **kwargs)


def post_process(func):
    ''' Used as a decorator to apply L{PostProcess} functions. '''
    def wrapper(*args, **kwargs):
        success = func(*args, **kwargs)
        if success:
            process_wrapper(*args, **kwargs)
        return success
    return wrapper


def pre_process(func):
    ''' Used as a decorator to apply L{PostProcess} functions. '''
    def wrapper(*args, **kwargs):
        process_wrapper(*args, **kwargs)
        success = func(*args, **kwargs)
        return success
    return wrapper


class PostProcess(object):

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

        def get_new_pmids(pmids, disease_code):
            ''' Find PMIDs in a list that are not in the elastic index. '''
            chunk_size = 800
            pmids_found = []
            time.sleep(5)
            for i in range(0, len(pmids), chunk_size):
                pmids_slice = pmids[i:i+chunk_size]
                query = ElasticQuery(Query.terms("PMID", pmids_slice), sources=['PMID', 'tags'])
                docs = Search(query, idx=section['index'], size=chunk_size).search().docs
                json_data = ''
                idx_name = ''
                idx_type = ''
                for doc in docs:
                    disease = getattr(doc, 'tags')['disease']
                    if disease is None:
                        disease = []
                    if disease_code not in disease:
                        # update disease attribute
                        disease.append(disease_code)
                        idx_name = doc._meta['_index']
                        idx_type = doc.type()

                        doc_data = {"update": {"_id": doc._meta['_id'], "_type": idx_type,
                                               "_index": idx_name, "_retry_on_conflict": 3}}
                        json_data += json.dumps(doc_data) + '\n'
                        json_data += json.dumps({'doc': {'tags': {'disease': disease}}}) + '\n'

                    pmids_found.append(getattr(doc, 'PMID'))

                if json_data != '':
                    Loader().bulk_load(idx_name, idx_type, json_data)
            return [pmid for pmid in pmids if pmid not in pmids_found]

        section_name = args[1]
        section_dir_name = args[2]
        base_dir_path = args[3]
        stage_dir = os.path.join(base_dir_path, 'STAGE', section_dir_name)

        if not os.path.exists(stage_dir):
            os.makedirs(stage_dir)

        section = kwargs['section']
        if 'output' not in section:
            return
        stage_file = os.path.join(stage_dir, section['output'] + '.json')
        download_file = os.path.join(base_dir_path, 'DOWNLOAD', section_dir_name, section['output'])

        tree = ET.parse(download_file)
        idlist = tree.find("IdList")
        ids = list(idlist.iter("Id"))
        pmids = [i.text for i in ids]

        parts = section_name.rsplit(':', 1)
        disease_code = parts[1].lower()

        if Search().index_exists(section['index']):
            pmids = get_new_pmids(pmids, disease_code)

        Pubs.fetch_details(pmids, stage_file, disease_code)


class IniParser(object):
    ''' Provides utility functions for ini file parsing. '''

    def read_ini(self, ini_file):
        ''' Download data defined in the ini file. '''
        # check for ini file in data pipeline home
        if not os.path.isfile(ini_file):
            DOWNLOAD_BASE_DIR = os.path.dirname(os.path.dirname(__file__))
            tmp = os.path.join(DOWNLOAD_BASE_DIR, 'data_pipeline', ini_file)
            if os.path.isfile(tmp):
                ini_file = tmp
        config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        config.read(ini_file)
        return config

    def process_sections(self, config, base_dir_path, sections=None):
        ''' Loop over all sections in the config file and process. '''
        success = False
        for section_name in config.sections():
            if sections is not None and section_name not in sections:
                continue

            section_dir_name = self._inherit_section(section_name, config)
            dir_path = os.path.join(base_dir_path, self.__class__.__name__.upper(), section_dir_name)
            success = self.process_section(section_name, section_dir_name, base_dir_path,
                                           dir_path=dir_path, section=config[section_name],
                                           stage=self.__class__.__name__)
        return success

    def process_section(self, section_name, section_dir_name, base_dir_path,
                        dir_path='.', section=None, stage=None):
        raise NotImplementedError("Inheriting class should implement this  method")

    def _inherit_section(self, section_name, config):
        ''' Add in parameters from another config section when a double colon
        is found in the name. '''
        if '::' in section_name:
            inherit = section_name.split('::', maxsplit=1)[0]
            if inherit in config:
                for key in config[inherit]:
                    config[section_name][key] = config[inherit][key]
            return inherit
        return section_name
