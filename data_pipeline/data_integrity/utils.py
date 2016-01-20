''' Utils for data_integrity '''
import logging
from elastic.query import Query, ScoreFunction, FunctionScoreQuery
from elastic.search import ElasticQuery, Search
import random
import requests
import json
from elastic.elastic_settings import ElasticSettings
logger = logging.getLogger(__name__)


class DataIntegrityUtils(object):

    @classmethod
    def get_rdm_feature_id(cls, idx, idx_type, qbool=Query.match_all(), sources=[], field=None):
        ''' Get a random feature id from the indices. DEPRECATED USE IN DJANGO-ELASTIC '''
        doc = cls.get_rdm_docs(idx, idx_type, qbool=qbool, sources=sources, size=1)[0]

        if field is not None:
            return getattr(doc, field)

        return doc.doc_id()

    @classmethod
    def get_rdm_docs(cls, idx, idx_type, qbool=Query.match_all(), sources=[], size=1):
        ''' Get a random doc from the indices. DEPRECATED USE IN DJANGO-ELASTIC '''
        score_function1 = ScoreFunction.create_score_function('random_score', seed=random.randint(0, 1000000))

        search_query = ElasticQuery(FunctionScoreQuery(qbool, [score_function1], boost_mode='replace'),
                                    sources=sources)
        elastic = Search(search_query=search_query, size=size, idx=idx, idx_type=idx_type)
        try:
            return elastic.search().docs
        except IndexError:
            return cls.get_rdm_docs(idx, idx_type, qbool, sources, size)

    @classmethod
    def get_rdm_feature_ids(cls, idx, idx_type, qbool=Query.match_all(), sources=[], field=None, size=1):
        ''' Get random feature_ids from the indices. DEPRECATED USE IN DJANGO-ELASTIC '''
        docs = cls.get_rdm_docs(idx, idx_type, qbool=qbool, sources=sources, size=size)

        ids = []
        for doc in docs:
            if field is not None:
                ids.append(getattr(doc, field))
            else:
                ids.append(doc.doc_id())

        return ids

    @classmethod
    def get_docs_count(cls, idx, idx_type):
        '''Get doc counts'''
        elastic = Search(idx=idx, idx_type=idx_type)
        return elastic.get_count()['count']

    @classmethod
    def fetch_from_ensembl(cls, gene_id):
        '''Lookup ensembl via restful call'''
        server = "http://rest.ensembl.org"
        ext = "/lookup/id/" + gene_id + "?content-type=application/json;expand=1;db_type=core;object_type=Gene"

        url = server+ext
        logger.debug(url)

        response = requests.get(url)

        if response.ok:
            data = json.loads(response.content.decode('utf-8'))
            return data
        else:
            return None

    @classmethod
    def fetch_xref_from_ensembl(cls, gene_id):
        '''Lookup ensembl entrez mapping via restful call'''
        server = "http://rest.ensembl.org"
        ext = "/xrefs/id/" + gene_id + "?content-type=application/json;expand=1;external_db=EntrezGene"

        url = server+ext
        logger.debug(url)

        response = requests.get(url)

        if response.ok:
            data = json.loads(response.content.decode('utf-8'))
            return data
        else:
            return None

    @classmethod
    def fetch_from_elastic(cls, idx, idx_type, feature_ids):
        '''Lookup pydgin elastic'''
        query = ElasticQuery(Query.ids(feature_ids))
        elastic = Search(query, idx=ElasticSettings.idx(idx, idx_type=idx_type), size=5)
        docs = elastic.search().docs
        return docs
