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
from elastic.query import Query, AndFilter, RangeQuery, BoolQuery
from elastic.elastic_settings import ElasticSettings

chr_band = '10p15.1'
tier_cutoff = 2
build = 38

idx = ElasticSettings.idx('REGION')
Delete.docs_by_query(idx, idx_type='disease_locus')

buildSort = {"sort": [
        {"build_info.start":
            {"order": "asc",
             "nested_path": "build_info",
             "nested_filter": {
                    "term": {
                        "build_info.build": build
                    }
                }
             }
         }]}


def add_disease_locus(seqid, locus_id, regionName, disease, tier, species, weight, doc_ids):
    data = {
        "region_name": disease+" "+regionName,
        "disease": disease,
        "tier": tier,
        "species": species,
        "tags": {"weight": weight},
        "locus_id": locus_id,
        "seqid": seqid,
        "hits": doc_ids
    }
    #    "suggest": {"input": [disease+" "+regionName, regionName], "weight": weight}
    resp = requests.put(ElasticSettings.url()+'/' + idx+'/disease_locus/'+locus_id, data=json.dumps(data))
    if resp.status_code != 201:
        print(str(resp.content))
        print("Problem loading "+getattr(doc, "disease")+" "+regionName)
    else:
        print("Loaded "+locus_id+" - "+regionName)

Search.index_refresh(idx)
diseases_by_seqid = Agg('diseases_by_seqid', 'terms', {"size": 0, "field": "disease"})
disease_hits = Agg('disease_hits', 'reverse_nested', {}, sub_agg=diseases_by_seqid)
seq_hits = Agg('seq_hits', 'terms', {'field': 'build_info.seqid', 'size': 0}, sub_agg=disease_hits)
build_info = Agg('build_info', 'nested', {"path": 'build_info'}, sub_agg=[seq_hits])

qnested = ElasticQuery(Query.nested('build_info', Query.term("build_info.build", build)))
elastic = Search(qnested, idx=idx, aggs=Aggs(build_info), search_type='count')
resultObj = elastic.search()
resultAggs = resultObj.aggs
seq_hits = getattr(resultAggs['build_info'], 'seq_hits')['buckets']

for chr_bucket in seq_hits:
    seqid = chr_bucket['key'].upper()
    
    for disease_bucket in chr_bucket['disease_hits']['diseases_by_seqid']['buckets']:
        # print(disease_bucket)
        disease_code = disease_bucket['key']
        print("\n"+str(seqid)+" - "+disease_code)
        query2 = ElasticQuery.filtered_bool(
            Query.nested("build_info", BoolQuery(
                must_arr=[Query.term("build_info.seqid", chr_bucket['key']), Query.term("build_info.build", build)])),
            BoolQuery(must_arr=[Query.term("disease", disease_code), RangeQuery("tier", lte=tier_cutoff)]),
            sources=["disease", "marker", "chr_band", "tier", "tags.weight", "species", "build_info"]
            )
        elastic2 = Search(search_query=query2, idx=idx, idx_type='hits',
                          size=int(disease_bucket['doc_count']+1), qsort=Sort(buildSort))
        results = elastic2.search()
        minPos = 0
        maxPos = 0
        tier = 4
        weight = 0
        regionCount = 1
        regionName = ''
        species = ''
        doc_ids = []
        if len(results.docs) > 0:                        
            for doc in results.docs:
                # print(doc)
                os.system("curl -XPOST '"+ElasticSettings.url()+"/"+idx+"/hits/" + doc.doc_id() +
                          "/_update?pretty' -d '{\"doc\": {\"disease_locus\": \"TBC\"}}' > /dev/null 2>&1")
                build_info = None
                for b in getattr(doc, 'build_info'):
                    if b['build'] == build:
                        build_info = b
                if build_info is None:
                    print("ERROR - no build information found for b"+str(build))
                    continue

                # print(getattr(doc, "disease")+"\t"+getattr(doc, "marker")+"\t" + getattr(doc, "chr_band") + "\t" +
                #      build_info['seqid'] + "\t" + str(build_info['start']) + "\t" + str(build_info['end']))
                if minPos == 0 and maxPos == 0:
                    minPos = build_info['start']
                    maxPos = build_info['end']

                if minPos <= build_info['start'] <= maxPos or minPos <= build_info['end'] <= maxPos:
                    minPos = min(minPos, build_info['start'])
                    maxPos = max(maxPos, build_info['end'])
                    tier = min(tier, int(getattr(doc, "tier")))
                    regionName = getattr(doc, "chr_band")
                    species = getattr(doc, "species")
                    weight = max(getattr(doc, "tags")['weight'], weight)
                    doc_ids.append(doc.doc_id())
                else:
                    locus_id = disease_code.upper()+"_"+str(seqid)+format(regionCount, '03d')
                    # print("id="+locus_id+"\tregion_name="+getattr(doc, "disease")+" "+regionName+"\ttier="+str(tier) +
                    #      "\tdisease="+getattr(doc, "disease")+"\tspecies="+species)
                    add_disease_locus(seqid, locus_id, regionName, disease_code.upper(), tier, species, weight, doc_ids)
                    for d in doc_ids:
                        os.system("curl -XPOST '"+ElasticSettings.url()+"/"+idx+"/hits/"+d+"/_update?pretty' -d '" +
                                  "{\"doc\": {\"disease_locus\": \""+locus_id+"\"}}' > /dev/null 2>&1")
                    regionCount += 1
                    doc_ids = []
                    doc_ids.append(doc.doc_id())
                    minPos = build_info['start']
                    maxPos = build_info['end']
                    tier = int(getattr(doc, "tier"))
                    regionName = getattr(doc, "chr_band")
                    species = getattr(doc, "species")
                    weight = getattr(doc, "tags")['weight']
            locus_id = disease_code.upper()+"_"+str(seqid)+format(regionCount, '03d')
            # print("id="+locus_id+"\tregion_name="+getattr(doc, "disease")+" "+regionName+"\ttier="+str(tier) +
            #      "\tdisease="+getattr(doc, "disease")+"\tspecies="+species)
            add_disease_locus(seqid, locus_id, regionName, disease_code.upper(), tier, species, weight, doc_ids)
            for d in doc_ids:
                os.system("curl -XPOST '"+ElasticSettings.url()+"/"+idx+"/hits/"+d+"/_update?pretty' -d '" +
                          "{\"doc\": {\"disease_locus\": \""+locus_id+"\"}}' > /dev/null 2>&1")