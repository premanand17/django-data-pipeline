import os
import sys
import requests
import json
import django

PYDGIN = None
if 'PYDGIN' in os.environ:
    PYDGIN = os.environ['PYDGIN']
else:
    print("ENV variable PYDGIN must be set to run this script")
    sys.exit()

sys.path.append(PYDGIN)
os.environ['DJANGO_SETTINGS_MODULE'] = 'pydgin.settings'
django.setup()

from elastic.aggs import Agg, Aggs
from elastic.search import ElasticQuery, FilteredQuery, Search, Sort, Delete, Update
from elastic.query import Query, Filter, AndFilter, RangeQuery, BoolQuery
from elastic.elastic_settings import ElasticSettings

chr_band = '10p15.1'
tier_cutoff = 2
build = 38

idx = ElasticSettings.idx('REGIONS')
Delete.docs_by_query(idx, idx_type='region')

def add_region(seqid, region_id, regionName, tier, species, weight, diseases, doc_ids):
    data = {
        "region_name": regionName,
        "tier": tier,
        "species": species,
        "tags": {"weight": weight, "disease": diseases},
        "region_id": region_id,
        "seqid": seqid,
        "disease_loci": doc_ids
    }

    resp = requests.put(ElasticSettings.url()+'/' + idx+'/region/'+region_id, data=json.dumps(data))
    if resp.status_code != 201:
        print(str(resp.content))
        print("Problem loading "+regionName)
    else:
        print("Loaded "+region_id+" - "+regionName)

Search.refresh(idx)
locus_start = Agg('locus_start', 'min', {'field': 'build_info.start'})
locus_end = Agg('locus_end', 'max', {'field': 'build_info.end'})
match_agg = Agg('filtered_result', 'filter', Query.match("build_info.build", build).query_wrap(),
                sub_agg=[locus_start, locus_end])
build_filter = Agg('build_filter', 'nested', {"path": 'build_info'},
                   sub_agg=[match_agg, locus_start])
locus_tier = Agg('locus_tier', 'min', {"field": "tier"})
locus_weight = Agg('locus_weight', 'max', {"field": "tags.weight"})
# elastic 1.7.X
# diseases_by_seqid = Agg('diseases_by_seqid', 'terms', {"size": 0, "field": "disease_locus",
#                                                       "order": {"build_filter>locus_start": "asc"}},
#                        sub_agg=[build_filter, locus_tier, locus_weight])
# elastic 2.x
diseases_by_seqid = Agg('diseases_by_seqid', 'terms', {"size": 0, "field": "disease_locus"},
                        sub_agg=[build_filter, locus_tier, locus_weight])
disease_hits = Agg('disease_hits', 'reverse_nested', {}, sub_agg=diseases_by_seqid)
seq_hits = Agg('seq_hits', 'terms', {'field': 'build_info.seqid', 'size': 0}, sub_agg=disease_hits)
build_info = Agg('build_info', 'nested', {"path": 'build_info'}, sub_agg=[seq_hits])

qnested = Query.nested('build_info', Query.term("build_info.build", build))
query = ElasticQuery.filtered(qnested, Filter(BoolQuery(must_not_arr=[Query.term("disease_locus", "TBC")])))
elastic = Search(query, idx=idx, aggs=Aggs(build_info), search_type='count')

resultObj = elastic.search()
resultAggs = resultObj.aggs
seq_hits = getattr(resultAggs['build_info'], 'seq_hits')['buckets']

''' loop round seqid buckets '''
for chr_bucket in seq_hits:
    seqid = chr_bucket['key'].upper()
    print("\n"+str(seqid))
    minPos = 0
    maxPos = 0
    tier = 4
    weight = 0
    regionCount = 1
    regionName = ''
    species = ''
    doc_ids = []
    diseases = []
    ''' loop round disease_locus buckets '''
    buckets = chr_bucket['disease_hits']['diseases_by_seqid']['buckets']
    sorted_buckets = sorted(buckets, key=lambda e: e['build_filter']['locus_start']['value'])
    # for locus in chr_bucket['disease_hits']['diseases_by_seqid']['buckets']:
    for locus in sorted_buckets:
        locus_id = locus['key']
        locus_start = int(locus['build_filter']['filtered_result']['locus_start']['value'])
        locus_end = int(locus['build_filter']['filtered_result']['locus_end']['value'])
        locus_tier = int(locus['locus_tier']['value'])
        locus_weight = int(locus['locus_weight']['value'])

        os.system("curl -XPOST '"+ElasticSettings.url()+"/"+idx+"/disease_locus/"+locus_id+"/_update?pretty' -d '" +
                  "{\"doc\": {\"region\": \"None\"}}' > /dev/null 2>&1")
        # Update.update_doc(part_doc='{"doc": {"region": "None"}}')

        if minPos == 0 and maxPos == 0:
            minPos = locus_start
            maxPos = locus_end

        if species == '':
            results = Search(search_query=ElasticQuery(Query.ids(locus_id), sources=['species']),
                             idx=idx, idx_type='disease_locus').search()
            species = getattr(results.docs[0], "species")

        if regionName == '':
            results = Search(search_query=ElasticQuery(Query.ids(locus_id), sources=['region_name']),
                             idx=idx, idx_type='disease_locus').search()
            regionName = getattr(results.docs[0], "region_name").split()[1]

        if minPos <= locus_start <= maxPos or minPos <= locus_end <= maxPos:
            minPos = min(minPos, locus_start)
            maxPos = max(maxPos, locus_end)
            tier = min(tier, locus_tier)
            weight = max(locus_weight, weight)
            doc_ids.append(locus_id)
        else:
            region_id = regionName+"_"+format(regionCount, '03d')
            results = Search(search_query=ElasticQuery(Query.ids(doc_ids), sources=['disease']),
                             idx=idx, idx_type='disease_locus').search()
            for d in results.docs:
                diseases.append(getattr(d, "disease"))
            add_region(seqid, region_id, regionName, tier, species, weight, diseases, doc_ids)
            for d in doc_ids:
                os.system("curl -XPOST '"+ElasticSettings.url()+"/"+idx+"/disease_locus/"+d+"/_update?pretty' -d '" +
                          "{\"doc\": {\"region\": \""+region_id+"\"}}' > /dev/null 2>&1")

            regionCount += 1

            results = Search(search_query=ElasticQuery(Query.ids(locus_id), sources=['species']),
                             idx=idx, idx_type='disease_locus').search()
            species = getattr(results.docs[0], "species")

            results = Search(search_query=ElasticQuery(Query.ids(locus_id), sources=['region_name']),
                             idx=idx, idx_type='disease_locus').search()
            regionName = getattr(results.docs[0], "region_name").split()[1]

            doc_ids = []
            diseases = []
            doc_ids.append(locus_id)
            minPos = locus_start
            maxPos = locus_end
            tier = locus_tier
            weight = locus_weight

    region_id = regionName+"_"+format(regionCount, '03d')
    results = Search(search_query=ElasticQuery(Query.ids(doc_ids), sources=['disease']),
                     idx=idx, idx_type='disease_locus').search()
    for d in results.docs:
        diseases.append(getattr(d, "disease"))
    add_region(seqid, region_id, regionName, tier, species, weight, diseases, doc_ids)
    for d in doc_ids:
        os.system("curl -XPOST '"+ElasticSettings.url()+"/"+idx+"/disease_locus/"+d+"/_update?pretty' -d '" +
                  "{\"doc\": {\"region\": \""+region_id+"\"}}' > /dev/null 2>&1")
