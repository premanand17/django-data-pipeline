=============
data_pipeline
=============

The data pipeline is configured using ini files. 

Publication Pipeline
--------------------

The publication.ini is used to define the configuration for downloading,
staging/processing and loading/indexing. To run the publication pipeline::

    ./manage.py publications --dir tmp --ini publications.ini \
                             --steps download stage load

Or step-by-step::

    ./manage.py publications --dir tmp --ini publications.ini \
                             --sections GENE --steps download stage
    ./manage.py publications --dir tmp --ini publications.ini \
                             --sections GENE --steps load

    ./manage.py publications --dir tmp --ini publications.ini \
                             --sections [DISEASE::T1D],[DISEASE::CRO] \
                             --steps  download load

The publication pipeline is incremental so that when run multiple times it
will query EUTILS only for new publications it finds and add those to the index.

Note a useful terms aggregation for finding the number of documents per disease::

    curl 'http://127.0.0.1:9200/publications/_search?size=1&from=0&pretty' \
       -d '{"aggs": {"disease_groups": {"terms": {"field": "disease", "size": 0}}}}'

    curl 'http://127.0.0.1:9200/publications/_search?size=0&from=0&pretty' \
       -d '{"aggs": {"missing_disease_groups": {"missing": {"field": "disease"}}}}'
