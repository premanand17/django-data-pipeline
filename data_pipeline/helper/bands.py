from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader, JSONLoader
import logging

logger = logging.getLogger(__name__)


class Bands(object):
    ''' Cytobands data '''

    @classmethod
    def mapping(cls, idx, idx_type):
        ''' Load the mapping for the cytobands index.
        seqid    - chromosome
        start    - start position
        stop     - stop position
        name     - cytoband name
        gieStain - Giemsa stain results. gneg, gpos50, gpos75, gpos25, gpos100, acen, gvar, stalk
        '''
        props = MappingProperties(idx_type)
        props.add_property("seqid", "string", index="not_analyzed") \
             .add_property("start", "integer") \
             .add_property("stop", "integer") \
             .add_property("name", "string", analyzer="full_name") \
             .add_property("giestain", "string", index="not_analyzed")

        ''' create index and add mapping '''
        load = Loader()
        options = {"indexName": idx, "shards": 2}
        load.mapping(props, idx_type, analyzer=Loader.KEYWORD_ANALYZER, **options)
        return props

    @classmethod
    def idx(cls, bands_f, idx, idx_type):
        ''' Parse and load data for cytobands. '''
        new_docs = []
        chunk_size = 450
        count = 0
        for band in bands_f:
            parts = band.strip().split('\t')
            doc = {}
            doc['seqid'] = parts[0].replace('chr', '')
            doc['start'] = int(parts[1])
            doc['stop'] = int(parts[2])
            doc['name'] = parts[3]
            doc['giestain'] = parts[4]
            doc['_id'] = doc['seqid'] + doc['name']
            new_docs.append(doc)
            if count > chunk_size:
                # check id's are current dbsnp rs id's
                JSONLoader().load(new_docs, idx, idx_type)
                new_docs = []
                count = 0
            count += 1
        if len(new_docs) > 0:
            JSONLoader().load(new_docs, idx, idx_type)


class Chrom(object):
    ''' Chromosome lengths. '''

    @classmethod
    def mapping(cls, idx, idx_type):
        ''' Load the mapping for the chromosome type in the bands index.
        seqid    - chromosome
        length   - sequence length
        '''
        props = MappingProperties(idx_type)
        props.add_property("seqid", "string", index="not_analyzed") \
             .add_property("length", "integer")

        ''' create index and add mapping '''
        load = Loader()
        options = {"indexName": idx, "shards": 1}
        load.mapping(props, idx_type, analyzer=Loader.KEYWORD_ANALYZER, **options)
        return props

    @classmethod
    def idx(cls, bands_f, idx, idx_type):
        ''' Parse and load data for chromosome lengths. '''
        new_docs = []
        chunk_size = 500
        count = 0
        for band in bands_f:
            parts = band.strip().split('\t')
            doc = {}
            doc['length'] = int(parts[1])
            doc['_id'] = parts[0].replace('chr', '')
            new_docs.append(doc)
            if count > chunk_size:
                # check id's are current dbsnp rs id's
                JSONLoader().load(new_docs, idx, idx_type)
                new_docs = []
                count = 0
            count += 1
        if len(new_docs) > 0:
            JSONLoader().load(new_docs, idx, idx_type)
