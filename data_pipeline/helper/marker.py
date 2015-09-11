from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader, JSONLoader
from elastic.query import TermsFilter, Query
from elastic.search import ElasticQuery, Search
from elastic.elastic_settings import ElasticSettings
import logging

logger = logging.getLogger(__name__)


class ImmunoChip(object):
    ''' Immunochip marker data '''

    @classmethod
    def ic_mapping(cls, idx, idx_type, test_mode=False):
        ''' Load the mapping for the immunochip marker index.
        id          - current rs marker id
        name        - IC alias
        allele_a/b  - alleles
        synonyms    - other synonyms
        build_info  - build number and position
        internal_id - internal marker alias
        strand      - strand
        is_par      - is in a pseudoautosomal region
        '''
        props = MappingProperties(idx_type)
        props.add_property("id", "string", analyzer="full_name") \
             .add_property("allele_a", "string", index="not_analyzed") \
             .add_property("allele_b", "string", index="not_analyzed") \
             .add_property("synonyms", "string", analyzer="full_name") \
             .add_property("build_info", "object") \
             .add_property("is_par", "string") \
             .add_property("name", "string", analyzer="full_name") \
             .add_property("internal_id", "integer") \
             .add_property("strand", "string", index="not_analyzed") \
             .add_property("suggest", "completion",
                           index_analyzer="full_name", search_analyzer="full_name")

        ''' create index and add mapping '''
        load = Loader()
        options = {"indexName": idx, "shards": 5}
        if not test_mode:
            load.mapping(props, idx_type, analyzer=Loader.KEYWORD_ANALYZER, **options)
        return props

    @classmethod
    def immunochip_mysql_2_idx(cls, ic_f, idx_name, idx_type):
        ''' Parse and load data for immunochip markers. '''
        new_docs = []
        chunk_size = 450
        count = 0
        current_marker_ids = []
        for ic in ic_f:
            parts = ic.strip().split('\t')
            if parts[0] == 'ilmn_id':
                continue

            doc = {}
            doc['allele_a'] = parts[1]
            doc['allele_b'] = parts[2]

            syns = set()
            syns.add(parts[0])
            current_marker_id = ''
            for m_id in parts[3:8]:
                if m_id != '\\N' and m_id != 'AMBIG':
                    current_marker_id = m_id
                    syns.add(m_id)
                else:
                    current_marker_id = ''
            if current_marker_id in syns:
                syns.remove(current_marker_id)

            suggests = []
            if len(syns) > 0:
                doc['synonyms'] = list(syns)
                suggests.extend(list(syns))
            if current_marker_id != '':
                doc['id'] = current_marker_id
                current_marker_ids.append(current_marker_id)
            doc['build_info'] = [{'build': '36', 'position': parts[8], 'seqid': parts[9]},
                                 {'build': '37', 'position': parts[10], 'seqid': parts[11]},
                                 {'build': '38', 'position': parts[12], 'seqid': parts[13]}]
            if parts[14] != '\\N':
                doc['is_par'] = parts[14]  # pseudoautosomal

            if parts[18] != '\\N':         # use marker_mart_141 if present
                doc['internal_id'] = int(parts[18])
            else:
                doc['internal_id'] = int(parts[15])
            doc['name'] = parts[16]
            if parts[17] != '\\N':
                doc['strand'] = parts[17]
            suggests.append(doc['name'])

            doc['suggest'] = {}
            doc['suggest']["input"] = suggests
            new_docs.append(doc)

            if count > chunk_size:
                # check id's are current dbsnp rs id's
                cls.check_rs_ids(current_marker_ids, new_docs)
                JSONLoader().load(new_docs, idx_name, idx_type)

                new_docs = []
                current_marker_ids = []
                count = 0
            count += 1
        if len(new_docs) > 0:
            JSONLoader().load(new_docs, idx_name, idx_type)

    @classmethod
    def check_rs_ids(cls, current_marker_ids, new_docs):
        ''' Check id's in immunochip docs are current in dbsnp. If not then
        move the id to the synonym fields. '''
        terms_filter = TermsFilter.get_terms_filter("id", current_marker_ids)
        query = ElasticQuery.filtered(Query.match_all(), terms_filter, sources='id')
        elastic = Search(query, idx=ElasticSettings.idx('MARKER', idx_type='MARKER'), size=len(current_marker_ids))
        marker_docs = elastic.search().docs

        not_current_marker_ids = []
        for m_id in current_marker_ids:
            if not cls._contains_id(marker_docs, m_id):
                not_current_marker_ids.append(m_id)
        if len(not_current_marker_ids) == 0:
            return

        ''' check rshigh if the marker id has merged, see docs:
        www.ncbi.nlm.nih.gov/projects/SNP/snp_db_table_description.cgi?t=RsMergeArch
        '''
        terms_filter = TermsFilter.get_terms_filter("rshigh", not_current_marker_ids)
        query = ElasticQuery.filtered(Query.match_all(), terms_filter, sources=['rscurrent', 'rshigh', "build_id"])
        elastic = Search(query, idx=ElasticSettings.idx('MARKER', idx_type='HISTORY'), size=len(current_marker_ids))
        history_docs = elastic.search().docs
        rshistory = {}
        for h_doc in history_docs:
            if getattr(h_doc, 'build_id') < 142:
                logger.error("MARKER MERGE BUILD < 142: " + getattr(h_doc, 'rshigh') + ' build: ' +
                             str(getattr(h_doc, 'build_id')))
            if getattr(h_doc, 'rscurrent') != 'rs':
                rshistory[getattr(h_doc, 'rshigh')] = getattr(h_doc, 'rscurrent')

        for m_id in not_current_marker_ids:
            # no longer a current marker id so move to synonym
            for n_doc in new_docs:
                if 'id' in n_doc and n_doc['id'] == m_id:
                    logger.debug("Marker no longer current: "+m_id)
                    if m_id in rshistory:
                        n_doc['id'] = rshistory[m_id]
                        n_doc['internal_id'] = int(rshistory[m_id].replace('rs', ''))
                    else:
                        del n_doc['id']
                    if m_id not in n_doc['synonyms']:
                        n_doc['synonyms'] = m_id

    @classmethod
    def _contains_id(cls, docs, mid):
        ''' Check for id in the elasticsearch docs. '''
        for doc in docs:
            if getattr(doc, 'id') == mid:
                return True
        return False
