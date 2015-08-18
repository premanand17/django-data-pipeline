''' Helper functions and classes '''
import os
import configparser
import time
import xml.etree.ElementTree as ET

from elastic.search import Search, ElasticQuery
from elastic.query import Query, TermsFilter
from .helper.pubs import Pubs
import json
from elastic.management.loaders.loader import Loader
import re
import gzip
import logging
from data_pipeline.helper.gene import Gene
from data_pipeline.helper.gene_interactions import GeneInteractions
from data_pipeline.helper.gene_pathways import GenePathways

# Get an instance of a logger
logger = logging.getLogger(__name__)


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
    def _get_stage_file(cls, *args, **kwargs):
        section = kwargs['section']
        section_dir_name = args[2]
        base_dir_path = args[3]
        stage_dir = os.path.join(base_dir_path, 'STAGE', section_dir_name)
        if not os.path.exists(stage_dir):
            os.makedirs(stage_dir)
        if 'output' in section:
            return os.path.join(stage_dir, section['output'] + '.json')
        elif 'files' in section:
            return os.path.join(stage_dir, section['files'] + '.json')

    @classmethod
    def _get_stage_output_file(cls, *args, **kwargs):
        section = kwargs['section']
        section_dir_name = args[2]
        base_dir_path = args[3]
        stage_dir = os.path.join(base_dir_path, 'STAGE', section_dir_name)
        if not os.path.exists(stage_dir):
            os.makedirs(stage_dir)
        if 'output' in section:
            return os.path.join(stage_dir, section['stage_output'])
        elif 'files' in section:
            return os.path.join(stage_dir, section['files'] + '.out')

    @classmethod
    def _get_download_file(cls, *args, **kwargs):
        section = kwargs['section']
        section_dir_name = args[2]
        base_dir_path = args[3]
        if 'output' in section:
            return os.path.join(base_dir_path, 'DOWNLOAD', section_dir_name, section['output'])
        elif 'files' in section:
            download_files = []
            ff = section['files'].split(',')
            if len(ff) > 1:
                for f in ff:
                    download_file = os.path.join(base_dir_path, 'DOWNLOAD', section_dir_name,
                                                 f.strip())
                    download_files.append(download_file)
                    print(download_files)
                return download_files
            else:
                return os.path.join(base_dir_path, 'DOWNLOAD', section_dir_name, section['files'])

    ''' Pipeline methods '''
    @classmethod
    def ensembl_gene_parse(cls, *args, **kwargs):
        ''' Parse gene GTF file from ensembl. '''
        stage_file = cls._get_stage_file(*args, **kwargs)
        download_file = cls._get_download_file(*args, **kwargs)
        Gene.gene_mapping(kwargs['section']['index'], kwargs['section']['index_type'])
        with gzip.open(download_file, 'rt') as ensembl_gene_f:
            with open(stage_file, 'w') as outfile:
                json.dump(Gene.ensembl_gene_parse(ensembl_gene_f), outfile, indent=0)

    @classmethod
    def ensmart_gene_parse(cls, *args, **kwargs):
        ''' Parse result from ensembl mart. '''
        download_file = cls._get_download_file(*args, **kwargs)
        with open(download_file, 'rt') as ensmart_f:
            Gene.ensmart_gene_parse(ensmart_f, kwargs['section']['index'])

    @classmethod
    def gene2ensembl_parse(cls, *args, **kwargs):
        ''' Parse gene2ensembl file from NCBI. '''
        download_file = cls._get_download_file(*args, **kwargs)
        with gzip.open(download_file, 'rt') as gene2ens_f:
            Gene.gene2ensembl_parse(gene2ens_f, kwargs['section']['index'])

    @classmethod
    def gene_info_parse(cls, *args, **kwargs):
        ''' Parse gene_info file from NCBI. '''
        download_file = cls._get_download_file(*args, **kwargs)
        idx = kwargs['section']['index']

        with gzip.open(download_file, 'rt') as gene_info_f:
            Gene.gene_info_parse(gene_info_f, idx)

    @classmethod
    def gene_pub_parse(cls, *args, **kwargs):
        ''' Parse gene2pubmed file from NCBI. '''
        download_file = cls._get_download_file(*args, **kwargs)
        with gzip.open(download_file, 'rt') as gene_pub_f:
            Gene.gene_pub_parse(gene_pub_f, kwargs['section']['index'])

    ''' Publication methods '''
    @classmethod
    def get_new_pmids(cls, pmids, idx, disease_code=None):
        ''' Find PMIDs in a list that are not in the elastic index. '''
        chunk_size = 800
        pmids_found = set()
        pmids_found_add = pmids_found.add
        time.sleep(5)

        for i in range(0, len(pmids), chunk_size):
            pmids_slice = pmids[i:i+chunk_size]
            terms_filter = TermsFilter.get_terms_filter("PMID", pmids_slice)
            query = ElasticQuery.filtered(Query.match_all(), terms_filter, sources=['PMID', 'tags'])

            docs = Search(query, idx=idx, size=chunk_size).search().docs
            json_data = ''

            for doc in docs:
                pmids_found_add(getattr(doc, 'PMID'))
                if disease_code is not None:
                    tags = getattr(doc, 'tags')
                    if 'disease' in tags:
                        disease = tags['disease']
                    else:
                        disease = []
                    if disease_code not in disease:
                        # update disease attribute
                        disease.append(disease_code)
                        tags['disease'] = disease
                        idx_name = doc._meta['_index']
                        idx_type = doc.type()

                        doc_data = {"update": {"_id": doc._meta['_id'], "_type": idx_type,
                                               "_index": idx_name, "_retry_on_conflict": 3}}
                        json_data += json.dumps(doc_data) + '\n'
                        json_data += json.dumps({'doc': {'tags': tags}}) + '\n'

            if json_data != '':
                Loader().bulk_load(idx_name, idx_type, json_data)

        return [pmid for pmid in pmids if pmid not in pmids_found]

    @classmethod
    def unique(cls, *args, **kwargs):
        ''' Combine a list of compressed files. '''
        section = kwargs['section']
        stage_file = cls._get_stage_file(*args, **kwargs)
        download_file = cls._get_download_file(*args, **kwargs)

        pmids = set()
        with gzip.open(download_file, 'rt') as outf:
            seen_add = pmids.add
            for x in outf:
                if not x.startswith('9606\t'):
                    continue
                pmid = re.split('\t', x)[2].strip()
                if pmid not in pmids:
                    seen_add(pmid)
        new_pmids = cls.get_new_pmids(list(pmids), section['index'])
        print(len(new_pmids))
        Pubs.fetch_details(new_pmids, stage_file)

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
    def gene_interaction_parse(cls, *args, **kwargs):
        stage_output_file = cls._get_stage_output_file(*args, **kwargs)
        download_file = cls._get_download_file(*args, **kwargs)
        section = kwargs['section']
        GeneInteractions.gene_interaction_parse(download_file, stage_output_file, section)

    @classmethod
    def gene_pathway_parse(cls, *args, **kwargs):
        stage_output_file = cls._get_stage_output_file(*args, **kwargs)
        download_files = cls._get_download_file(*args, **kwargs)
        section = kwargs['section']
        GenePathways.gene_pathway_parse(download_files, stage_output_file, section)

    @classmethod
    def xmlparse(cls, *args, **kwargs):
        section_name = args[1]
        section = kwargs['section']
        stage_file = cls._get_stage_file(*args, **kwargs)
        download_file = cls._get_download_file(*args, **kwargs)
        if download_file is None:
            return

        tree = ET.parse(download_file)
        idlist = tree.find("IdList")
        ids = list(idlist.iter("Id"))
        pmids = [i.text for i in ids]
        npmids = len(pmids)

        parts = section_name.rsplit(':', 1)
        disease_code = parts[1].lower()

        if Search().index_exists(section['index']):
            pmids = cls.get_new_pmids(pmids, section['index'], disease_code=disease_code)

        logger.debug("Total No. of PMIDs in "+args[1]+": "+str(npmids))
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
