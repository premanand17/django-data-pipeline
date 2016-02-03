from elastic.management.loaders.loader import DelimeterLoader, MappingProperties
from elastic.management.loaders.loader import Loader


class RecombinationRates():
    ''' HapMap Recombination Rates data '''

    @classmethod
    def idx(cls, f, idx, idx_type):
        ''' Parse and load data for HapMap Recombination Data. '''
        column_names = ["seqid", "position", "recombination_rate", "genetic_map_position"]
        loader = DelimeterLoader()
        loader.mapping_json = None
        loader.load(column_names, f, idx, idx_type)

    @classmethod
    def mapping(cls, idx, idx_type):
        ''' Load the mapping for the recombination rates index.
        seqid    - chromosome
        position
        recombination_rate
        genetic_map_position
        '''
        props = MappingProperties(idx_type)
        props.add_property("seqid", "string", index="not_analyzed") \
             .add_property("position", "integer", index="not_analyzed") \
             .add_property("recombination_rate", "float", index="not_analyzed") \
             .add_property("genetic_map_position", "float", index="not_analyzed")

        ''' create index and add mapping '''
        load = Loader()
        options = {"indexName": idx, "shards": 2}
        load.mapping(props, idx_type, **options)
        return props
