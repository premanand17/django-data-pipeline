''' Test for LD data. '''
from django.test import TestCase
from django.core.urlresolvers import reverse
from elastic.query import Query, BoolQuery, RangeQuery
from data_pipeline.data_integrity.utils import DataIntegrityUtils
from elastic.elastic_settings import ElasticSettings
from random import randint
import json
import requests
import random
import logging

logger = logging.getLogger(__name__)


class LDTest(TestCase):

    def test_ld_list_population(self):
        ''' Test data from the restful LD server against the Ensembl restful
        LD interface (http://rest.ensembl.org/documentation/info/ld_id_get). '''
        m1 = self._get_random_marker()
        populations = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']
        pop = populations[randint(0, len(populations)-1)]

        url = reverse('rest:ld-list') + "?format=json&m1=%s&dataset=%s&window_size=25000&rsq=0.8" % (m1, pop)
        response = self.client.get(url, format='json')
        json_content = json.loads(response.content.decode("utf-8"))[0]

        if 'error' in json_content and json_content['error'] is not None:
            logger.warn(m1+' :: '+json_content['error'])
            return
        lds = json_content['ld']

        server = "http://rest.ensembl.org"
        ext = "/ld/human/%s?population_name=1000GENOMES:phase_3:%s&r2=0.8" % (m1, pop)
        ens_lds = requests.get(server+ext, headers={"Content-Type": "application/json"}).json()
        if 'error' in ens_lds:
            logger.warn(m1+' :: '+ens_lds['error'])
            return

        for ld in lds:
            m2 = ld['marker2']
            for ens_ld in ens_lds:
                if 'variation2' in ens_ld and m2 == ens_ld['variation2']:
                    r2_1 = round(float(ens_ld['r2']), 2)
                    r2_2 = ld['rsquared']
                    dprime_1 = round(float(ens_ld['d_prime']), 2)
                    dprime_2 = round(float(ld['dprime']), 2)

                    logger.debug(pop+' :: '+m1+' '+m2+' '+str(r2_1)+' '+str(r2_2)+' '+str(dprime_1)+' '+str(dprime_2))
                    self.assertEqual(r2_1, r2_2, m2+' r2:'+str(r2_1)+' '+str(r2_2))
                    self.assertEqual(dprime_1, dprime_2, m2+' dprime1:'+str(dprime_1)+' '+str(dprime_2))

        if len(lds) - len(ens_lds) != 0:
            logger.warn("DIFFERENCE IN LD FOR "+m1)

    def _get_random_marker(self):
        (idx, idx_type) = ElasticSettings.idx('MARKER', 'MARKER').split('/')

        seqid = random.randint(1, 10)
        qbool = BoolQuery(must_arr=[Query.term("seqid", seqid), RangeQuery("tags.weight", gte=80)])
        doc = DataIntegrityUtils.get_rdm_docs(idx, idx_type, qbool=qbool, sources=['id', 'start'], size=1)[0]
        return getattr(doc, 'id')
