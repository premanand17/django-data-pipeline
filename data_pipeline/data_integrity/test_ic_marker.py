''' Data integrity tests. '''
from django.test import TestCase
from elastic.elastic_settings import ElasticSettings
from elastic.search import ScanAndScroll, ElasticQuery, Search
from elastic.result import Document
import logging
from elastic.query import Query, TermsFilter
import requests
import sys

logger = logging.getLogger(__name__)


class ImmunoChipMarkerDataTest(TestCase):
    '''IC marker test '''
    neg_internal_id = 0

    def _rs_exists(self, rsid):
        url = "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=" + rsid.replace("rs", "")
        r = requests.get(url)
        if r.status_code != 200:
            return False
        if 'invalid snp_id' in r.text:
            return False
        return True

    def test_positions(self):
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
                            msg = ("DIFFERENT POSITIONS ID: "+str(doc_internal_id)+":\t" +
                                   str(getattr(doc1, "name"))+": "+pos1['position']+" ("+doc1.doc_id()+")\t" +
                                   str(getattr(doc2, "name"))+": "+pos2['position']+" ("+doc2.doc_id()+")\t")
                            try:
                                terms_filter = TermsFilter.get_terms_filter("start", [pos1['position'],
                                                                                      pos2['position']])
                                query = ElasticQuery.filtered(Query.term("seqid", pos1['seqid']), terms_filter)
                                elastic = Search(query, idx=ElasticSettings.idx('MARKER', 'MARKER'))
                                docs_by_pos = elastic.search().docs
                                found = False
                                for d in docs_by_pos:
                                    msg += getattr(d, "id")+": "+str(getattr(d, "start"))+"\t"
                                    if getattr(d, "id") == 'rs'+str(doc_internal_id):
                                        found = True

                                if not found:
                                    msg += 'rs'+str(doc_internal_id)
                                    if self._rs_exists('rs'+str(doc_internal_id)):
                                        msg += ' EXISTS IN DBSNP\t'
                                    else:
                                        msg += ' NOT IN DBSNP\t'
                                logger.error(msg)
                            except KeyError:
                                logger.error(msg)
                    internal_id[doc_internal_id].append(doc1)
                else:
                    internal_id[doc_internal_id] = [doc1]

        ScanAndScroll.scan_and_scroll(ElasticSettings.idx('MARKER', idx_type='IC'), call_fun=check_hits)
        print("LEN = "+str(len(internal_id)))

    def test_positions2(self):
        ''' Test the rs position matches dbsnp position for SNVs. '''
        def check_hits(resp_json):
            rsids = {}
            docs = [Document(hit) for hit in resp_json['hits']['hits']]
            for doc in docs:
                rsid = getattr(doc, "id")
                if rsid is not None:
                    rsids[rsid] = doc
            rsids_keys = list(rsids.keys())
            terms_filter = TermsFilter.get_terms_filter("id", rsids_keys)
            query = ElasticQuery.filtered(Query.match_all(), terms_filter)
            elastic = Search(query, idx=ElasticSettings.idx('MARKER', 'MARKER'), size=len(rsids_keys))
            docs_by_rsid = elastic.search().docs
            for doc in docs_by_rsid:
                info = getattr(doc, "info")
                if 'VC=SNV' not in info:
                    continue
                rsid = getattr(doc, "id")
                pos1 = getattr(doc, "start")
                ic_doc = rsids[rsid]
                pos2 = self._get_highest_build(ic_doc)['position']
                if abs(int(pos1) - int(pos2)) > 1:
                    is_par = getattr(ic_doc, 'is_par')
                    allele_a = getattr(ic_doc, 'allele_a')
                    if is_par is None and not (allele_a == 'D' or allele_a == 'I'):
                        logger.error("CHECK IC POSITIONS: "+getattr(ic_doc, 'name') +
                                     ' '+str(pos2)+" "+rsid+' '+str(pos1))

        ScanAndScroll.scan_and_scroll(ElasticSettings.idx('MARKER', idx_type='IC'), call_fun=check_hits)

    def test_internal_ids(self):
        ''' Test that the internal id matches the rs id. '''
        def check_hits(resp_json):
            docs = [Document(hit) for hit in resp_json['hits']['hits']]
            for doc1 in docs:
                internal_id = getattr(doc1, "internal_id")
                rsid = getattr(doc1, "id")
                if (internal_id is None and rsid < 0) or rsid is None:
                    continue
                rsid = rsid.replace('rs', '')
                self.assertEqual(int(internal_id), int(rsid), str(rsid)+" ::: "+str(internal_id))

        ScanAndScroll.scan_and_scroll(ElasticSettings.idx('MARKER', idx_type='IC'), call_fun=check_hits)

    def test_neg_internal_ids(self):
        ''' Test that the internal id matches the rs id. '''
        def check_hits(resp_json):
            docs = [Document(hit) for hit in resp_json['hits']['hits']]
            for doc1 in docs:
                internal_id = getattr(doc1, "internal_id")
                if int(internal_id) < self.neg_internal_id:
                    self.neg_internal_id = int(internal_id)

        ScanAndScroll.scan_and_scroll(ElasticSettings.idx('MARKER', idx_type='IC'), call_fun=check_hits)
        print(self.neg_internal_id)

    def _get_highest_build(self, doc):
        builds = getattr(doc, "build_info")
        high_build = {'build': '0'}
        for build in builds:
            if int(high_build['build']) < int(build['build']):
                high_build = build
        return high_build
