''' Utils for data_integrity '''
import logging
from elastic.query import Query
from elastic.search import ElasticQuery, Search
import requests
import json
from elastic.elastic_settings import ElasticSettings
logger = logging.getLogger(__name__)


class DataIntegrityUtils(object):

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
