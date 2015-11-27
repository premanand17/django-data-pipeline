''' Data integrity tests for publications index '''
from django.test import TestCase
from data_pipeline.utils import IniParser
from data_pipeline.download import HTTPDownload, FTPDownload
import xml.etree.ElementTree as ET
from urllib.parse import urljoin
from elastic.query import Query, BoolQuery, Filter
from elastic.search import ElasticQuery, Search
from elastic.elastic_settings import ElasticSettings
import data_pipeline
import logging
import os
import gzip
import re

logger = logging.getLogger(__name__)


class DiseasePublicationTest(TestCase):
    ''' Disease publication tests. '''
    TEST_DATA_DIR = os.path.dirname(data_pipeline.__file__) + '/tests/data'
    DISEASES = []

    @classmethod
    def setUpClass(cls):
        ''' Retrieve the publication list for each disease from NCBI. '''
        ini = IniParser()
        config = ini.read_ini('publications.ini')
        res = Search(ElasticQuery(Query.match_all(), sources=['code']), idx=ElasticSettings.idx('DISEASE')).search()
        sections = ''
        for doc in res.docs:
            sections += 'DISEASE::'+getattr(doc, 'code').upper()+','
        # sections = 'DISEASE::T1D,DISEASE::MS,DISEASE::SLE'

        # download ncbi publication lists for each disease
        for section_name in config.sections():
            if sections is not None and not ini._is_section_match(section_name, sections):
                continue
            ini._inherit_section(section_name, config)
            logger.debug(section_name)
            section = config[section_name]
            disease = section_name.split('::')[1]
            file_name = 'disease_pub_'+disease+'.tmp'
            HTTPDownload().download(section['location']+"?"+section['http_params'],
                                    cls.TEST_DATA_DIR, file_name=file_name)
            DiseasePublicationTest.DISEASES.append(disease)
        print()

    @classmethod
    def tearDownClass(cls):
        ''' Remove disease publication files. '''
        super(DiseasePublicationTest, cls)
        for disease in cls.DISEASES:
            file_name = 'disease_pub_'+disease+'.tmp'
            os.remove(os.path.join(cls.TEST_DATA_DIR, file_name))

    def test_pub_disease_counts(self):
        ''' Check all publications exist in the publication index. '''
        for disease in DiseasePublicationTest.DISEASES:
            pmids = self._get_pmids(disease)
            disease_code = disease.lower()
            elastic = Search(search_query=ElasticQuery(BoolQuery(b_filter=Filter(Query.ids(pmids)))),
                             idx=ElasticSettings.idx('PUBLICATION'), size=len(pmids)*2)
            self.assertEquals(elastic.get_count()['count'], len(pmids), 'Count for '+disease_code)

            # check for differences in pmids
            docs = elastic.search().docs
            pmids_in_idx = [getattr(doc, 'pmid') for doc in docs]
            pmids_diff = list(set(pmids) - set(pmids_in_idx))
            self.assertEquals(len(pmids_diff), 0)

    def test_pubs_disease_tags(self):
        ''' Check the number of disease publications against the number of tags.disease and
        report differences`. '''
        count = True
        msg = ''
        for disease in DiseasePublicationTest.DISEASES:
            pmids = self._get_pmids(disease)
            disease_code = disease.lower()
            elastic = Search(search_query=ElasticQuery(BoolQuery(
                         b_filter=Filter(Query.term('tags.disease', disease_code))), sources=['pmid']),
                         idx=ElasticSettings.idx('PUBLICATION'), size=len(pmids)*2)
            res = elastic.get_count()
            msg += disease_code+'\tINDEX: '+str(res['count'])+'\tNCBI: '+str(len(pmids))
            if res['count'] != len(pmids):
                count = False
                docs = elastic.search().docs
                pmids_in_idx = [getattr(doc, 'pmid') for doc in docs]
                pmids_diff1 = [pmid for pmid in pmids_in_idx if pmid not in pmids]
                pmids_diff2 = [pmid for pmid in pmids if pmid not in pmids_in_idx]
                if len(pmids_diff1) > 0:
                    msg += '\textra PMIDs: '+str(pmids_diff1)
                if len(pmids_diff2) > 0:
                    msg += '\tmissing PMIDs: '+str(pmids_diff2)
            msg += '\n'

        print(msg)
        self.assertTrue(count, 'Count for disease tags')

    def _get_pmids(self, disease):
        ''' Get the PMID list from NCBI XML file. '''
        xmlfile = 'disease_pub_'+disease+'.tmp'
        tree = ET.parse(os.path.join(DiseasePublicationTest.TEST_DATA_DIR, xmlfile))
        idlist = tree.find("IdList")
        ids = list(idlist.iter("Id"))
        return [i.text for i in ids]


class GenePublicationTest(TestCase):
    '''Gene publication tests. '''

    NUM_DIFF = 1000  # difference allowed for between index and NCBI

    def test_gene_pubs(self):
        ''' Check the difference between the pubs indexed and those from the gene_pub file
        from the NCBI. If the publication pipeline has not been run recently there is likely
        to be a difference. This is allowed for with the NUM_DIFF variable. If there is a
        larger difference than this then the publication pipeline should be run. '''
        ini = IniParser()
        config = ini.read_ini('publications.ini')
        section = config['GENE']

        file_name = 'gene_pub_test.tmp'
        download_file = os.path.join(DiseasePublicationTest.TEST_DATA_DIR, file_name)
        success = FTPDownload().download(urljoin(section['location'], section['files']),
                                         DiseasePublicationTest.TEST_DATA_DIR, file_name=file_name)
        self.assertTrue(success, 'downloaded gene publications file')

        pmids = set()
        with gzip.open(download_file, 'rt') as outf:
            seen_add = pmids.add
            for x in outf:
                if not x.startswith('9606\t'):
                    continue
                pmid = re.split('\t', x)[2].strip()
                if pmid not in pmids:
                    seen_add(pmid)

        elastic = Search(search_query=ElasticQuery(BoolQuery(b_filter=Filter(Query.ids(list(pmids)))),
                                                   sources=['pmid']),
                         idx=ElasticSettings.idx('PUBLICATION'), size=len(pmids)*2)
        self.assertLess(len(pmids)-elastic.get_count()['count'], GenePublicationTest.NUM_DIFF,
                        'Count for gene publications')

        # check for differences in pmids
        docs = elastic.search().docs
        pmids_in_idx = [getattr(doc, 'pmid') for doc in docs]
        pmids_diff = list(pmids - set(pmids_in_idx))
        self.assertLess(len(pmids_diff), GenePublicationTest.NUM_DIFF)
        os.remove(download_file)
