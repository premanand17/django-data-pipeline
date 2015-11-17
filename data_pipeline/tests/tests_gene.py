''' Tests for the staging. '''
from django.test import TestCase
from django.core.management import call_command
from django.utils.six import StringIO
import os
import data_pipeline
import shutil
import logging
from data_pipeline.helper.gene_interactions import GeneInteractions
import json
import re
from data_pipeline.helper.gene import Gene
from data_pipeline.utils import IniParser
from data_pipeline.helper.gene_pathways import GenePathways
logger = logging.getLogger(__name__)


def tearDownModule():
    '''This is module level tearDown, will run after all the test are ran in this module
    Removes the test data dirs
    '''
    logging.debug('Module teardown...')
    app_data_dir = os.path.dirname(data_pipeline.__file__)
    test_data_dir = app_data_dir + '/tests/data'
    stage_data_dir = test_data_dir + '/STAGE'
    if(os.path.exists(stage_data_dir)):
        shutil.rmtree(stage_data_dir)


class GeneInteractionStagingTest(TestCase):
    '''Test interaction staging'''

    def setUp(self):
        '''Runs before each of the tests run from this class..creates the tests/data dir'''
        logging.debug('Class setup...' + self.__class__.__name__)
        self.out = StringIO()
        self.ini_file = os.path.join(os.path.dirname(__file__), 'test_download.ini')
        self.app_data_dir = os.path.dirname(data_pipeline.__file__)
        self.test_data_dir = self.app_data_dir + '/tests/data'
        call_command('pipeline', '--steps', 'stage', sections='INTACT',
                     dir=self.test_data_dir, ini=self.ini_file, stdout=self.out)
        logging.debug(self.test_data_dir)

    def test_intact(self):
        ''' Test intact staging. '''
        logging.debug("Test intact staging")

        self.assertTrue(os.path.exists(self.test_data_dir + '/STAGE/INTACT'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/INTACT/intact.zip.out'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/INTACT/intact.zip.json'))
        self.assertTrue(os.access(self.test_data_dir + '/STAGE/INTACT/intact.zip.json', os.R_OK))

        # check if the .out file is in right format - this file will have interactors with ensemblid
        with open(self.test_data_dir + '/STAGE/INTACT/intact.zip.out', 'r') as f:
            header = f.readline()
            first_line = f.readline()
            f.close()

        self.assertEqual(header, "interactorA\tinteractorB\tpubmed\n", "header is ok")
        self.assertEqual(first_line, "ENSG00000078053\tENSG00000159082\t10542231\n", "first line is ok")

        # check if the json outputput file has the right data
        with open(self.test_data_dir + '/STAGE/INTACT/intact.zip.json', 'r') as f:
            docs_line = f.readline()
            interactor_line = f.readline()
            f.close()

        self.assertIn("docs", docs_line, "docs line OK")
        self.assertIn("interaction_source", interactor_line, "contains interaction_source")
        self.assertIn("intact", interactor_line, "contains intact as interaction_source")

    def test_intact_staged_output_files(self):
        '''Test to check the processed json output files'''
        # preinitialize dict
        expected_interactors = [
            {"interactor": "ENSG00000115641", "pubmed": "11001931"},
            {"interactor": "ENSG00000136068", "pubmed": "9437013"},
            {"interactor": "ENSG00000196924", "pubmed": "9437013"},
            {"interactor": "ENSG00000171552", "pubmed": "10446169"},
            {"interactor": "ENSG00000075142", "pubmed": "10748169"},
            {"interactor": "ENSG00000135018", "pubmed": "11076969"},
            {"interactor": "ENSG00000185043", "pubmed": "10366599"}
            ]

        parent = r'"_parent": "ENSG00000143801"'
        matched_json_str = None
        with open(self.test_data_dir + '/STAGE/INTACT/intact.zip.json', 'r') as f:
            for line in f:
                if re.search(parent, line):
                    matched_json_str = line.rstrip()
                    matched_json_str = matched_json_str[:-1]
                    break

        json_data = json.loads(matched_json_str)
        parent = json_data["_parent"]
        stage_interactors = json_data["interactors"]
        interaction_source = json_data["interaction_source"]

        self.assertEqual(parent, "ENSG00000143801", "Parent is ENSG00000143801")
        self.assertEqual(interaction_source, "intact", "Source is intact")

        expected_dict = dict((i['interactor'], i['pubmed']) for i in expected_interactors)
        staged_dict = dict((i['interactor'], i['pubmed']) for i in stage_interactors)
        # assertDictEqual should return True even if the keys-value pairs are in different order
        self.assertDictEqual(expected_dict, staged_dict, "Staging OK")

    def test_bioplex(self):
        ''' Test bioplex staging. '''
        logging.debug("Test bioplex staging")
        call_command('pipeline', '--steps', 'stage', sections='BIOPLEX',
                     dir=self.test_data_dir, ini=self.ini_file, stdout=self.out)

        self.assertTrue(os.path.exists(self.test_data_dir + '/STAGE/BIOPLEX'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.out'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.json'))
        self.assertTrue(os.access(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.json', os.R_OK))

        # check if the .out file is in right format - this file will have interactors with ensemblid
        with open(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.out', 'r') as f:
            header = f.readline()
            first_line = f.readline()
            f.close()

        self.assertEqual(header, "interactorA\tinteractorB\n", "header is ok")
        self.assertEqual(first_line, "ENSG00000196839\tENSG00000196604\n", "first line is ok")

        # check if the json outputput file has the right data
        with open(self.test_data_dir + '/STAGE/BIOPLEX/BioPlex_interactionList_v4.tsv.json', 'r') as f:
            docs_line = f.readline()
            interactor_line = f.readline()
            f.close()

        self.assertIn("docs", docs_line, "docs line OK")
        self.assertIn("interaction_source", interactor_line, "contains interaction_source")
        self.assertIn("bioplex", interactor_line, "contains bioplex as interaction_source")


class GenePathwayStagingTest(TestCase):
    '''Test gene pathway staging'''
    def setUp(self):
        '''Runs before each of the tests run from this class..creates the tests/data dir'''
        logging.debug('Class setup...' + self.__class__.__name__)
        self.out = StringIO()
        self.ini_file = os.path.join(os.path.dirname(__file__), 'test_download.ini')
        self.app_data_dir = os.path.dirname(data_pipeline.__file__)
        self.test_data_dir = self.app_data_dir + '/tests/data'
        logging.debug(self.test_data_dir)

    def test_msigdb(self):
        ''' Test msigdb staging. '''
        logging.debug("Test msigdb staging")
        call_command('pipeline', '--steps', 'stage', sections='MSIGDB',
                     dir=self.test_data_dir, ini=self.ini_file, stdout=self.out)

        self.assertTrue(os.path.exists(self.test_data_dir + '/STAGE/MSIGDB'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.biocarta.v5.0.entrez.gmt.json'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.kegg.v5.0.entrez.gmt.json'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.reactome.v5.0.entrez.gmt.json'))
        self.assertTrue(os.path.isfile(self.test_data_dir + '/STAGE/MSIGDB/c5.all.v5.0.entrez.gmt.json'))
        self.assertTrue(os.access(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.biocarta.v5.0.entrez.gmt.json', os.R_OK))

        pathway = r'"REACTOME_APOPTOTIC_CLEAVAGE_OF_CELLULAR_PROTEINS"'
        matched_json_str = None
        with open(self.test_data_dir + '/STAGE/MSIGDB/c2.cp.reactome.v5.0.entrez.gmt.json', 'r') as f:
            for line in f:
                if re.search(pathway, line):
                    matched_json_str = line.rstrip()
                    matched_json_str = matched_json_str[:-1]
                    break

        json_data = json.loads(matched_json_str)
        self.assertEqual(json_data['source'], 'reactome', 'Got reactome as source')
        self.assertEqual(json_data['pathway_name'], 'REACTOME_APOPTOTIC_CLEAVAGE_OF_CELLULAR_PROTEINS',
                         'Got right pathway name')
        self.assertEquals(len(json_data['gene_sets']), 38, "Found 38 Genes in gene_sets")


class GeneInteractionProcessTest(TestCase):
    '''Test functions in GeneInteractions class'''

    def test_check_binary_interactions(self):
        '''Test if the grouping of the binary interactions including the evidence is working correctly'''
        gene_interactors = {
            'gene0': [{'gene1:evidence1'}, {'gene2:evidence2'}],
            'gene1': [{'gene0:evidence0'}, {'gene2:evidence2'}, {'gene3:evidence3'}]
            }

        interactorA = "gene0"
        interactorB = "gene4"
        (gene_interactors_list_a, gene_interactors_list_b) = GeneInteractions._check_binary_interactions(gene_interactors,  # @IgnorePep8
                                                                                                         interactorA,
                                                                                                         interactorB,
                                                                                                         evidence_id="evidence4") # @IgnorePep8

        gene_interactors[interactorA] = gene_interactors_list_a
        gene_interactors[interactorB] = gene_interactors_list_b

        processed_interactors1 = {'gene0': [{'gene1:evidence1'}, {'gene2:evidence2'}, {'gene4': 'evidence4'}],
                                  'gene1': [{'gene0:evidence0'}, {'gene2:evidence2'}, {'gene3:evidence3'}],
                                  'gene4': [{'gene0': 'evidence4'}],
                                  }

        self.assertDictEqual(gene_interactors, processed_interactors1, "Grouping OK first iteration")

        interactorA = "gene4"
        interactorB = "gene5"
        (gene_interactors_list_a, gene_interactors_list_b) = GeneInteractions._check_binary_interactions(gene_interactors, # @IgnorePep8
                                                                                                         interactorA,
                                                                                                         interactorB,
                                                                                                         evidence_id="evidence5") # @IgnorePep8

        gene_interactors[interactorA] = gene_interactors_list_a
        gene_interactors[interactorB] = gene_interactors_list_b

        processed_interactors2 = {'gene0': [{'gene1:evidence1'}, {'gene2:evidence2'}, {'gene4': 'evidence4'}],
                                  'gene1': [{'gene0:evidence0'}, {'gene2:evidence2'}, {'gene3:evidence3'}],
                                  'gene4': [{'gene0': 'evidence4'}, {'gene5': 'evidence5'}],
                                  'gene5': [{'gene4': 'evidence5'}],
                                  }

        self.assertDictEqual(gene_interactors, processed_interactors2, "Grouping OK second iteration")

    def test_group_binary_interactions(self):
        '''
        Test if the interactors are grouped correctly
        '''
        binary_interactions = [('0', '1'), ('0', '2'), ('1', '2'), ('2', '1'), ('1', '3'), ('2', '3'), ('3', '4'),
                               ('4', '5'), ('5', '6'), ('5', '7'), ('6', '8'), ('7', '8'), ('8', '9'), ('8', '9')]

        expected_dict = {
            '0': [{'1': None}, {'2': None}],
            '1': [{'0': None}, {'2': None}, {'3': None}],
            '2': [{'0': None}, {'1': None}, {'3': None}],
            '3': [{'1': None}, {'2': None}, {'4': None}],
            '4': [{'3': None}, {'5': None}],
            '5': [{'4': None}, {'6': None}, {'7': None}],
            '6': [{'5': None}, {'8': None}],
            '7': [{'5': None}, {'8': None}],
            '8': [{'6': None}, {'7': None}, {'9': None}],
            '9': [{'8': None}]
            }

        grouped_interactions = GeneInteractions._group_binary_interactions(binary_interactions)

        self.assertDictEqual(grouped_interactions, expected_dict, "Grouping of interactions done OK")
        self.assertTrue(10 == len(grouped_interactions), "dict contains 10 groups")
        self.assertIn('8', expected_dict, "key 8 is in dict")
        self.assertIn({'9': None}, grouped_interactions['8'], "value 9 is in dict['8'] ")
        self.assertIn({'1': None}, grouped_interactions['0'], "value 1 is in dict['0'] ")

    def test_interactor_json_decorator(self):
        '''Test if the dict passed is returned as json with the interactor and evidence as the keys'''
        json_a = GeneInteractions.interactor_json_decorator({"geneA": None})
        self.assertDictEqual(json_a, {'interactor': 'geneA'}, "JSON equal for geneA")

        json_b = GeneInteractions.interactor_json_decorator({"geneA": "evidenceA"})
        self.assertDictEqual(json_b, {'interactor': 'geneA', 'pubmed': 'evidenceA'},
                             "JSON equal for geneA with evidence")

    def test_interaction_json_decorator(self):
        '''Test if the interaction with interactors and source is returned correctly'''
        json_interaction = GeneInteractions.interaction_json_decorator("intact",
                                                                       "geneA",
                                                                       [{"geneB": "1234"},
                                                                        {"geneC": "2345"},
                                                                        {"geneD": "6789"}])

        self.assertIn('"_parent": "geneA"', json_interaction, "_parent found in json")
        self.assertIn('"interaction_source": "intact"', json_interaction, "_parent found in json")

        json_interaction = GeneInteractions.interaction_json_decorator("bioplex",
                                                                       "geneA",
                                                                       [{"geneB": None},
                                                                        {"geneC": None},
                                                                        {"geneD": None}])
        self.assertJSONEqual(json_interaction,
                             json.dumps({"interactors": [{"interactor": "geneB"},
                                                         {"interactor": "geneC"},
                                                         {"interactor": "geneD"}],
                                         "interaction_source": "bioplex", "_parent": "geneA"}),
                             "JSON for interaction equal")


class GeneConversionTest(TestCase):

    def test__check_gene_history(self):
        '''Test if the right newid is fetched from genehistory'''
        config = IniParser().read_ini("tests/test_download.ini")
        section = config["BIOPLEX"]
        self.assertIsNotNone(section, "Section is not none")

        gene_sets = ['339457', '197215', '26191']
        (newgene_ids, discontinued_ids) = Gene._check_gene_history(gene_sets, section)
        self.assertTrue(len(newgene_ids) == 1, "Got back one new id")
        self.assertIn('339457', newgene_ids, "Got back 339457 in new gene ids")
        self.assertTrue(len(discontinued_ids) == 1, "Got back one discontinued geneid")

    def test__replace_oldids_with_newids(self):
        '''Test if the old ids are getting replaced wih newids'''
        gene_sets = ['339457', '197215', '26191']
        new_gene_ids = {'339457': '85452'}
        replaced_gene_sets = Gene._replace_oldids_with_newids(gene_sets, new_gene_ids)
        self.assertEqual(replaced_gene_sets, ['85452', '197215', '26191'], "Replaced 339457 with 85452")

        discontinued_ids = ['197215']
        replaced_gene_sets = Gene._replace_oldids_with_newids(gene_sets, new_gene_ids, discontinued_ids)
        print(replaced_gene_sets)
        self.assertEqual(replaced_gene_sets, ['85452', '26191'], "Replaced 339457 with 85452")


class GenePathwayProcessTest(TestCase):

    def setUp(self):
        '''Runs before each of the tests run from this class..creates the tests/data dir'''
        self.ini_file = os.path.join(os.path.dirname(__file__), 'test_download.ini')
        self.app_data_dir = os.path.dirname(data_pipeline.__file__)
        self.test_data_dir = self.app_data_dir + '/tests/data'

        self.config = IniParser().read_ini("tests/test_download.ini")
        self.section = self.config["MSIGDB"]
        self.assertIsNotNone(self.section, "Section is not none")

    def test__get_pathway_source(self):
        '''Test if the right source name is returned given a download file '''
        download_file_kegg = '/dunwich/scratch/prem/tmp/download/DOWNLOAD/MSIGDB/c2.cp.kegg.v5.0.entrez.gmt'
        source = GenePathways._get_pathway_source(download_file_kegg)
        self.assertTrue(source == "kegg", "Got back kegg as source")

        download_file_reactome = '/dunwich/scratch/prem/tmp/download/DOWNLOAD/MSIGDB/c2.cp.reactome.v5.0.entrez.gmt'
        source = GenePathways._get_pathway_source(download_file_reactome)
        self.assertTrue(source == "reactome", "Got back reactome as source")

        download_file_biocarta = '/dunwich/scratch/prem/tmp/download/DOWNLOAD/MSIGDB/c2.cp.biocarta.v5.0.entrez.gmt'
        source = GenePathways._get_pathway_source(download_file_biocarta)
        self.assertTrue(source == "biocarta", "Got back biocarta as source")

        download_file_go = '/dunwich/scratch/prem/tmp/download/DOWNLOAD/MSIGDB/c5.all.v5.0.entrez.gmt'
        source = GenePathways._get_pathway_source(download_file_go)
        self.assertTrue(source == "GO", "Got back go as source")
