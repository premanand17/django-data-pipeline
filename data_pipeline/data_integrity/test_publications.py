''' Data integrity tests for publications index '''
from django.test import TestCase
from data_pipeline.utils import IniParser
from data_pipeline.download import HTTPDownload
import xml.etree.ElementTree as ET
from elastic.query import Query, BoolQuery, Filter
from elastic.search import ElasticQuery, Search
from elastic.elastic_settings import ElasticSettings
import data_pipeline
import logging
import os

logger = logging.getLogger(__name__)


class PublicationTest(TestCase):
    ''' Publication tests. '''

    def test_publication_disease_counts(self):
        ''' Retrieve the publication list for each disease from NCBI and check
        they all exist in the publication index.
        '''
        ini = IniParser()
        config = ini.read_ini('publications.ini')
        res = Search(ElasticQuery(Query.match_all(), sources=['code']), idx=ElasticSettings.idx('DISEASE')).search()
        sections = ''
        for doc in res.docs:
            sections += 'DISEASE::'+getattr(doc, 'code').upper()+','
        print(sections)

        TEST_DATA_DIR = os.path.dirname(data_pipeline.__file__) + '/tests/data'

        for section_name in config.sections():
            if sections is not None and not ini._is_section_match(section_name, sections):
                continue
            ini._inherit_section(section_name, config)
            logger.debug(section_name)
            section = config[section_name]

            file_name = 'disease_pub_'+section_name.split('::')[1]+'.tmp'
            self.assertTrue(HTTPDownload().download(section['location']+"?"+section['http_params'],
                                                    TEST_DATA_DIR, file_name=file_name))

            tree = ET.parse(os.path.join(TEST_DATA_DIR, file_name))
            idlist = tree.find("IdList")
            ids = list(idlist.iter("Id"))
            pmids = [i.text for i in ids]
            parts = section_name.rsplit(':', 1)
            disease_code = parts[1].lower()
            res = Search(search_query=ElasticQuery(BoolQuery(b_filter=Filter(Query.ids(pmids)))),
                         idx=ElasticSettings.idx('PUBLICATION')).get_count()
            self.assertEquals(res['count'], len(pmids), 'Count for '+disease_code)
            os.remove(os.path.join(TEST_DATA_DIR, file_name))
