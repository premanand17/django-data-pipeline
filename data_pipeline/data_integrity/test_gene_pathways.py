''' Data integrity tests for gene pathway index '''
from django.test import TestCase
import logging
from data_pipeline.utils import IniParser
import data_pipeline
import os
from data_pipeline.download import HTTPDownload
from elastic.elastic_settings import ElasticSettings
from data_pipeline.data_integrity.utils import DataIntegrityUtils
from elastic.query import BoolQuery, Query
from data_pipeline.helper.gene_pathways import GenePathways
import re
from data_pipeline.helper.gene import Gene
import shutil

logger = logging.getLogger(__name__)

TEST_DATA_DIR = os.path.dirname(data_pipeline.__file__) + '/tests/data/MSIGDB_TEST'


def setUpModule():
        config = IniParser().read_ini("download.ini")

        if not os.path.exists(TEST_DATA_DIR):
            os.makedirs(TEST_DATA_DIR)

        # Download MSIGDB file from source
        section = config["MSIGDB"]
        username = section['username'] if 'username' in section else None
        password = section['password'] if 'password' in section else None

        if 'files' in section:
                files = section['files'].split(",")
                for f in files:
                    f = f.replace(' ', '')
                    success = HTTPDownload.download(section['location']+"/"+f.strip(), TEST_DATA_DIR,
                                                    file_name=f, username=username, password=password)
                    logger.warn("Download status for msigdb " + str(success) + TEST_DATA_DIR)


def tearDownModule():
    # rm test dir
    shutil.rmtree(TEST_DATA_DIR)


class GenePathwayDataTest(TestCase):
    '''Gene Pathway Data Integrity test '''

    def test_gene_pathways(self):
        idx_key = 'GENE'
        idx_type_key = 'PATHWAY'

        idx = ElasticSettings.idx(idx_key, idx_type_key)
        (idx, idx_type) = idx.split('/')

        # Test doc count
        doc_count = DataIntegrityUtils.get_docs_count(idx, idx_type)
        self.assertGreater(doc_count, 2500, 'Gene doc count greater than 2500')

        # Get pathway doc - passing the pathway source (kegg, reactome, biocarta, go) and id . Also test with random id
        pathway_doc_kegg = self.get_pathway_doc("kegg")
        logger.debug(pathway_doc_kegg.__dict__)
        self.check_kegg_data(pathway_doc_kegg)

    def check_kegg_data(self, pathway_doc):
        # get the attributes from elastic doc
        pw_name_es = getattr(pathway_doc, 'pathway_name')
        pw_url_es = getattr(pathway_doc, 'pathway_url')
        gene_sets_es = set(getattr(pathway_doc, 'gene_sets'))

        config = IniParser().read_ini("download.ini")
        section = config["MSIGDB"]
        download_files = section['files'].split(",")

        # get the latest source file
        source_file = None
        for file_ in download_files:
            file_ = file_.replace(' ', '')
            source = GenePathways._get_pathway_source(file_)
            print(source)
            if source == 'kegg':
                source_file = file_
                break

        my_regex = r"^" + re.escape(pw_name_es) + r"\b"

        match_found = False
        pw_name_kegg = None
        pw_url_kegg = None
        pw_gene_sets_kegg = []
        if os.path.isfile(TEST_DATA_DIR + '/' + source_file):
            with open(TEST_DATA_DIR + '/' + source_file, "r") as data:
                for line in data:
                    if re.search(my_regex, line):
                        match_found = True
                        tmp_list = line.split()
                        pw_name_kegg = tmp_list[0]
                        pw_url_kegg = tmp_list[1]
                        pw_gene_sets_kegg = tmp_list[2:]

            if not match_found:
                self.assertTrue(1 == 2, 'No matching pathway found in file')
        else:
            logger.critical('File doesnt exists' + source_file)
            self.assertTrue(1 == 2, 'File doesnt exist ' + source_file)

        # compare the name and url
        self.assertEqual(pw_name_es, pw_name_kegg, 'pathway name in index and source same')
        self.assertEqual(pw_url_es, pw_url_kegg, 'pathway url in index and source same')

        # For genesets do a entrez id lookup
        section = config["ENSEMBL_GENE"]
        entrez_ensembl_dict = Gene._entrez_ensembl_lookup(pw_gene_sets_kegg, section, config)

        ensembl_list_kegg = [ensembl_id for entrez_id, ensembl_id in entrez_ensembl_dict.items()]  # @UnusedVariable
        ensembl_list_kegg = set(ensembl_list_kegg)

        diff = set()
        if(len(gene_sets_es) == len(ensembl_list_kegg)):
            # check if list contents are equal
            import collections
            self.assertTrue(collections.Counter(gene_sets_es) == collections.Counter(ensembl_list_kegg),
                            "two sets have equal elements")
        else:
            # find the missing one - Subtract.
            logger.warning('len of gene_sets_es ' + str(len(gene_sets_es)) +
                           ' len of gene_sets_pydgin ' + str(len(ensembl_list_kegg)))
            if(len(gene_sets_es) > len(ensembl_list_kegg)):
                diff = set(gene_sets_es) - set(ensembl_list_kegg)
            else:
                diff = set(ensembl_list_kegg) - set(gene_sets_es)

            # now check if these ids exists in history
            # Do a entrez to ensembl id lookup in gene history
            (newgene_ids, discontinued_ids) = Gene._check_gene_history(list(diff), config)  # @UnusedVariable
            logger.debug(newgene_ids)
            logger.debug(discontinued_ids)
            if(len(diff) == len(discontinued_ids)):
                self.assertEqual(len(diff), len(discontinued_ids),
                                 "The missing ids where found in gene_history as discontinued ids")

    def get_pathway_doc(self, pathway_source='kegg'):
        '''Fetch random and specific genes from elastic'''
        idx_key = 'GENE'
        idx_type_key = 'PATHWAY'

        idx = ElasticSettings.idx(idx_key, idx_type_key)
        (idx, idx_type) = idx.split('/')
        qbool_pw = BoolQuery().should([Query.term("source", pathway_source)])

        # Get random doc or specific if id is passed in query
        random_doc = DataIntegrityUtils.get_rdm_docs(idx, idx_type, qbool=qbool_pw, sources=[], size=1)
        doc = random_doc[0]

        return doc
