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
            for doc1 in docs:
                doc_internal_id = getattr(doc1, "internal_id")
                if doc_internal_id in internal_id:
                    pos1 = self._get_highest_build(doc1)
                    for doc2 in internal_id[doc_internal_id]:
                        pos2 = self._get_highest_build(doc2)
                        if pos2['position'] != pos1['position']:
                            msg = ("ID "+str(doc_internal_id)+" has different positions:\t" +
                                   str(getattr(doc1, "name"))+": "+pos1['position']+" ("+doc1.doc_id()+")\t" +
                                   str(getattr(doc2, "name"))+": "+pos2['position']+" ("+doc2.doc_id()+")\t")
                            try:
                                terms_filter = TermsFilter.get_terms_filter("start", [pos1['position'],
                                                                                      pos2['position']])
                                query = ElasticQuery.filtered(Query.term("seqid", pos1['seqid']), terms_filter)
                                elastic = Search(query, idx=ElasticSettings.idx('MARKER', 'MARKER'))
                                docs_by_pos = elastic.search().docs
                                for d in docs_by_pos:
                                    msg += getattr(d, "id")+": "+str(getattr(d, "start"))+"\t"
                                logger.error(msg)
                            except KeyError:
                                logger.error(msg)
                    internal_id[doc_internal_id].append(doc1)
                else:
                    internal_id[doc_internal_id] = [doc1]

        ScanAndScroll.scan_and_scroll(ElasticSettings.idx('MARKER', idx_type='IC'), call_fun=check_hits)
        print("LEN = "+str(len(internal_id)))

    def _get_highest_build(self, doc):
        builds = getattr(doc, "build_info")
        high_build = {'build': '0'}
        for build in builds:
            if int(high_build['build']) < int(build['build']):
                high_build = build
        return high_build
