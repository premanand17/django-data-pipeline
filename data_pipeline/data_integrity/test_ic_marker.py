''' Data integrity tests. '''
from django.test import TestCase
from elastic.elastic_settings import ElasticSettings
from elastic.search import ScanAndScroll, ElasticQuery, Search
from elastic.result import Document
import logging
from elastic.query import Query, TermsFilter

logger = logging.getLogger(__name__)


class ImmunoChipMarkerDataTest(TestCase):
    '''IC marker test '''

    def test_internal_ids(self):
        internal_id = {}

        def check_hits(resp_json):
            self.assertTrue('hits' in resp_json, 'scan and scroll hits')
            self.assertGreaterEqual(len(resp_json['hits']['hits']), 1)
            docs = [Document(hit) for hit in resp_json['hits']['hits']]
            for doc in docs:
                doc_internal_id = getattr(doc, "internal_id")
                if doc_internal_id in internal_id:
                    position1 = self._get_highest_build_pos(doc)
                    for doc2 in internal_id[doc_internal_id]:
                        position2 = self._get_highest_build_pos(doc2)
                        if position2 != position1:
                            try:
                                terms_filter = TermsFilter.get_terms_filter("id", getattr(doc, "synonyms"))
                                query = ElasticQuery.filtered(Query.match_all(), terms_filter)
                                elastic = Search(query, idx=ElasticSettings.idx('MARKER', 'MARKER'))
                                docs = elastic.search().docs
                                if len(docs) == 1:
                                    rs_current_id = getattr(docs[0], "id")
                                    rs_position = getattr(docs[0], "start")

                                    logger.error("ID "+str(doc_internal_id)+" has different positions:\t" +
                                                 str(getattr(doc, "name"))+": "+position1+"\t" +
                                                 str(getattr(doc2, "name"))+": "+position2+"\t" +
                                                 rs_current_id+": "+str(rs_position))
                                else:
                                    logger.error("ID "+str(doc_internal_id)+" has different positions:\t" +
                                                 str(getattr(doc, "name"))+": "+position1+"\t" +
                                                 str(getattr(doc2, "name"))+": "+position2)
                            except KeyError:
                                pass

                    internal_id[doc_internal_id].append(doc)
                else:
                    internal_id[doc_internal_id] = [doc]

        ScanAndScroll.scan_and_scroll(ElasticSettings.idx('MARKER', idx_type='IC'), call_fun=check_hits)
        print("LEN = "+str(len(internal_id)))

    def _get_highest_build_pos(self, doc):
        builds = getattr(doc, "build_info")
        high_build = {'build': '0'}
        for build in builds:
            if int(high_build['build']) < int(build['build']):
                high_build = build
        return high_build['position']
