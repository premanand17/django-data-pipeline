''' Data integrity tests for gene interaction index '''
from django.test import TestCase
from elastic.elastic_settings import ElasticSettings
import logging
from elastic.query import Query, BoolQuery
from data_pipeline.data_integrity.utils import DataIntegrityUtils
from data_pipeline.helper.gene import Gene
from data_pipeline.utils import IniParser
from data_pipeline.download import HTTPDownload, FTPDownload
import re
import zipfile
import csv
import data_pipeline
import os
logger = logging.getLogger(__name__)

TEST_DATA_DIR = os.path.dirname(data_pipeline.__file__) + '/tests/data'
INTACT = 'intact.zip'
BIOPLEX = 'bioplex.tmp'


def setUpModule():
        config = IniParser().read_ini("download.ini")
        # Download intact file from source
        section_intact = config["INTACT"]
        file_url = section_intact['location'] + section_intact['files']
        status = FTPDownload.download(file_url, TEST_DATA_DIR, INTACT)
        logger.warn("Download status for intact " + str(status))

        section_bioplex = config["BIOPLEX"]
        file_url = section_bioplex['location'] + section_bioplex['files']

        status = HTTPDownload.download(file_url, TEST_DATA_DIR, BIOPLEX)
        logger.warn("Download status for bioplex " + str(status))


def tearDownModule():
    os.remove(os.path.join(TEST_DATA_DIR, INTACT))
    os.remove(os.path.join(TEST_DATA_DIR, BIOPLEX))


class GeneInteractionDataTest(TestCase):
    '''Gene Interaction test '''

    def test_gene_interactions(self):
        '''Fetch random genes from elastic and compare the same with the results fetched directly from intact'''
        # elastic doc example:
        # "_source":{"interaction_source": "intact", "interactors": [
        # {"interactor": "ENSG00000206053", "pubmed": "16169070"},
        # {"interactor": "ENSG00000101474", "pubmed": "16169070"},
        # {"interactor": "ENSG00000065361", "pubmed": "16169070"},
        # {"interactor": "ENSG00000085465", "pubmed": "16169070"}]}

        idx_key = 'GENE'
        idx_type_key = 'INTERACTIONS'

        idx = ElasticSettings.idx(idx_key, idx_type_key)
        (idx, idx_type) = idx.split('/')

        # Test doc count
        doc_count = DataIntegrityUtils.get_docs_count(idx, idx_type)
        self.assertGreater(doc_count, 23000, 'Gene doc count greater than 60000')

        # Get interaction doc - passing the interaction source and id . Also test with random id
        (child_doc_bioplex, parent_doc_bioplex) = self.get_interaction_doc("bioplex", parent_id="ENSG00000164933")
        self.check_bioplex_data(child_doc_bioplex, parent_doc_bioplex)

        (child_doc_bioplex, parent_doc_bioplex) = self.get_interaction_doc("bioplex")
        self.check_bioplex_data(child_doc_bioplex, parent_doc_bioplex)

        (child_doc_intact, parent_doc_intact) = self.get_interaction_doc("intact", parent_id="ENSG00000188786")
        self.check_intact_data(child_doc_intact, parent_doc_intact)

        (child_doc_intact, parent_doc_intact) = self.get_interaction_doc("intact")
        self.check_intact_data(child_doc_intact, parent_doc_intact)

    def get_interaction_doc(self, interaction_source='intact', parent_id=None):
        '''Fetch random and specific genes from elastic'''
        idx_key = 'GENE'
        idx_type_key = 'INTERACTIONS'
        parent_idx_key = 'GENE'

        idx = ElasticSettings.idx(idx_key, idx_type_key)
        (idx, idx_type) = idx.split('/')

        if parent_id:
            qbool_intact = BoolQuery().must([Query.term("interaction_source", interaction_source),
                                            Query.term("_parent", parent_id)])
        else:
            qbool_intact = BoolQuery().should([Query.term("interaction_source", interaction_source)])

        # Get random doc or specific if id is passed in query
        docs_by_geneid = DataIntegrityUtils.get_rdm_docs(idx, idx_type, qbool=qbool_intact, sources=[], size=1)
        doc = docs_by_geneid[0]

        # Get parent doc
        parent_id = doc.parent()
        parent_docs = DataIntegrityUtils.fetch_from_elastic(idx_key, parent_idx_key, [parent_id])

        if parent_docs:
            self.assertTrue(len(parent_docs) >= 1, "Found 1 parent")
            parent_doc = parent_docs[0]
            return doc, parent_doc
        else:
            return self.get_interaction_doc("intact", parent_id)

    def check_intact_data(self, child_doc, parent_doc):
        '''Fetch interactors stored in elastic and compare them with what is fetched from source'''
        self.assertEqual(getattr(child_doc, "interaction_source"), 'intact', 'interaction_source is intact')

        # Get interactors already stored in our pipeline
        interactors = getattr(child_doc, 'interactors')
        pydgin_interactors = [interactor['interactor'] for interactor in interactors]

        parent_id = parent_doc.doc_id()
        self.assertEqual(parent_id, child_doc.parent(), 'Parent id ok')

        # Download intact file from source and search for the parent entrez id interactors
        if os.path.isfile(TEST_DATA_DIR + '/' + INTACT):
            parent_intact = set()
            zf = zipfile.ZipFile(TEST_DATA_DIR + '/' + INTACT, 'r')

            my_regex = re.escape(parent_id)
            if 'intact.txt' in zf.namelist():
                target_path = zf.extract(member='intact.txt', path=TEST_DATA_DIR)
                with open(target_path, encoding='utf-8') as csvfile:
                    reader = csv.reader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
                    for row in reader:
                        interactorA = row[22]
                        interactorB = row[23]

                        # line = '\t'.join(row)
                        matchA = re.search(my_regex, interactorA)
                        if matchA:
                            parent_intact.add(interactorB)

                        matchB = re.search(my_regex, interactorB)
                        if matchB:
                            parent_intact.add(interactorA)

                intact_interactors = set()
                for interactor_line in parent_intact:
                    match = re.search(r"(ENSG[0-9]*)", interactor_line)
                    if match:
                        result = match.group(0)
                        intact_interactors.add(result)

                self.assertEqual(len(pydgin_interactors), len(intact_interactors), "Interactors size equal")

                pydgin = list(pydgin_interactors)
                intact = list(intact_interactors)
                pydgin.sort()
                intact.sort()

                logger.debug(pydgin)
                logger.debug(intact)

                for x, y in zip(pydgin, intact):
                    self.assertEqual(x, y, 'Interactors are equal ' + x + '  <=> ' + y)

    def check_bioplex_data(self, child_doc, parent_doc):
        '''
        Get all interactors, collect the ensembl ids, convert them to entrez ids
        Fetch the source file from bioplex and search for the parent entrez id
        Compare if the interactors if count is same between two sets
        If there is difference, check if the entrez id is in gene_history
        '''
        config = IniParser().read_ini("download.ini")

        self.assertEqual(getattr(child_doc, "interaction_source"), 'bioplex', 'interaction_source is bioplex')

        # Get interactors
        interactors = getattr(child_doc, 'interactors')
        # Get ensembl ids
        ensembl_ids_interactors = [interactor['interactor'] for interactor in interactors]

        # Do a ensembl to entrez id lookup
        section = config["ENSEMBL_GENE"]
        ensembl_entrez_dict = Gene._ensembl_entrez_lookup(ensembl_ids_interactors, section)

        entrez_list_pydgin = set()
        for ensembl_id, entrez_id in ensembl_entrez_dict.items():  # @UnusedVariable
            entrez_list_pydgin.add(entrez_id)

        number_of_interactors_pydgin = len(interactors)

        parent_id = parent_doc.doc_id()
        self.assertEqual(parent_id, child_doc.parent(), 'Parent id ok')

        parent_entrez = getattr(parent_doc, "dbxrefs")["entrez"]
        # Download bioplex file from source and search for the parent entrez id interactors
        my_regex = r"\b" + re.escape(parent_entrez) + r"\b"
        interactor_counter = 0
        if os.path.isfile(TEST_DATA_DIR + '/' + BIOPLEX):
            entrez_list_bioplex = set()
            with open(TEST_DATA_DIR + '/' + BIOPLEX, "r") as data:
                for line in data:
                    if re.search(my_regex, line):
                        tmp_list = line.split()
                        if tmp_list[0] != parent_entrez:
                            entrez_list_bioplex.add(tmp_list[0])
                        if tmp_list[1] != parent_entrez:
                            entrez_list_bioplex.add(tmp_list[1])
                        interactor_counter += 1

        if(len(entrez_list_pydgin) == len(entrez_list_bioplex)):
            self.assertEqual(number_of_interactors_pydgin, interactor_counter,
                             "Interactor count is correct " + str(number_of_interactors_pydgin))
        else:
            # find the missing one - Subtract.
            diff = set()
            if(len(entrez_list_pydgin) > len(entrez_list_bioplex)):
                diff = entrez_list_pydgin - entrez_list_bioplex
            else:
                diff = entrez_list_bioplex - entrez_list_pydgin

            # now check if these ids exists in history
            # Do a entrez to ensembl id lookup in gene history
            (newgene_ids, discontinued_ids) = Gene._check_gene_history(list(diff), config)  # @UnusedVariable
            if(len(diff) == len(discontinued_ids)):
                self.assertEqual(len(diff), len(discontinued_ids),
                                 "The missing ids where found in gene_history as discontinued ids")
            else:
                # check in ensembl xref
                for gene_id in diff:
                    ensembl_gene_xref = DataIntegrityUtils.fetch_xref_from_ensembl(gene_id)
                    self.assertTrue(len(ensembl_gene_xref) == 0, "mapping not found in ensembl for " + gene_id)
