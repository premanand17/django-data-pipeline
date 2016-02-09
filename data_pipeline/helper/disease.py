''' Disease index data. '''
import json
import logging
import re

import requests

from elastic.elastic_settings import ElasticSettings
from elastic.management.loaders.loader import Loader
from elastic.management.loaders.mapping import MappingProperties


logger = logging.getLogger(__name__)


class Disease(object):
    ''' Disease data '''

    @classmethod
    def mapping(cls, idx, idx_type):
        ''' Create the mapping for disease indexing '''
        props = MappingProperties("disease")
        props.add_property("name", "string") \
             .add_property("code", "string") \
             .add_property("description", "string", index="not_analyzed") \
             .add_property("colour", "string", index="not_analyzed") \
             .add_property("tier", "integer", index="not_analyzed") \
             .add_property("suggest", "completion", analyzer="full_name")

        tags = MappingProperties("tags")
        tags.add_property("weight", "integer", index="not_analyzed")
        props.add_properties(tags)
        load = Loader()
        options = {"indexName": idx, "shards": 1}
        load.mapping(props, 'disease', analyzer=Loader.KEYWORD_ANALYZER, **options)

    @classmethod
    def idx(cls, disease_f, idx, idx_type):
        ''' Parse and load data for cytobands. '''
        for line in disease_f:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = re.split('\t', line)
            data = {
                "name": parts[0],
                "code": parts[2].lower(),
                "description": parts[1],
                "colour": parts[3],
                "tier": int(parts[4])
            }
            data['suggest'] = {}
            data['suggest']["input"] = [parts[2].lower(), parts[0]]
            data['suggest']["weight"] = 250
            resp = requests.put(ElasticSettings.url()+'/' +
                                idx+'/'+idx_type+'/'+parts[2].lower(),
                                data=json.dumps(data))
            if resp.status_code == 201:
                logger.debug("Loaded "+parts[0])
            else:
                logger.error("Problem loading "+parts[0])
