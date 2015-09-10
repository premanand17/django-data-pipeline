from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader, JSONLoader
from elastic.query import TermsFilter, Query
from elastic.search import ElasticQuery, Search
from elastic.elastic_settings import ElasticSettings
import sys


class ImmunoChip(object):
    ''' Immunochip marker data '''

    @classmethod
    def ic_mapping(cls, idx, idx_type, test_mode=False):
        ''' Load the mapping for the immunochip marker index. '''
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
        docs = []
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
            for m_id in parts[3:7]:
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
                doc['is_par'] = parts[14]
            doc['internal_id'] = int(parts[15])
            doc['name'] = parts[16]
            if parts[17] != '\\N':
                doc['strand'] = parts[17]
            suggests.append(doc['name'])

            doc['suggest'] = {}
            doc['suggest']["input"] = suggests
            docs.append(doc)

            if count > chunk_size:
#                 terms_filter = TermsFilter.get_terms_filter("id", current_marker_ids)
#                 query = ElasticQuery.filtered(Query.match_all(), terms_filter)
#                 elastic = Search(query, idx=ElasticSettings.idx('MARKER', 'immunochip'), size=chunk_size)
#                 marker_docs = elastic.search().docs
# 
#                 for m_id in current_marker_ids:
#                     found = False
#                     for m_doc in marker_docs:
#                         if getattr(m_doc, 'id') == m_id:
#                             found = True
#                     if not found:
#                         print(m_id)
#                 sys.exit()

                JSONLoader().load(docs, idx_name, idx_type)
                docs = []
                count = 0
            count += 1
        if len(docs) > 0:
            JSONLoader().load(docs, idx_name, idx_type)
