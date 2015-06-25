=============
data_pipeline
=============

The data pipeline is configured using ini files. 

Publication Pipeline
--------------------

The publication.ini is used to define the configuration for downloading,
staging/processing and loading/indexing. To run the publication pipeline:

    ./manage.py publications --dir tmp --ini publications.ini --steps  download stage  load

Or step-by-step:

    ./manage.py publications --dir tmp --ini publications.ini \
                             --sections GENE --steps download stage
    ./manage.py publications --dir tmp --ini publications.ini \
                             --sections GENE --steps load

    ./manage.py publications --dir tmp --ini publications.ini \
                             --sections [DISEASE::T1D],[DISEASE::CRO] --steps  download load

    
Note a useful terms aggregation for finding the number of documents per disease:

    curl 'http://127.0.01:9200/publications/_search?size=1&from=0&pretty' \
       -d '{"aggregations": {"test": {"terms": {"field": "disease", "size": 0}}}}'
