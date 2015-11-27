import logging
import sys
import csv
import re
import zipfile
import os
import json
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
from data_pipeline.helper.gene import Gene

logger = logging.getLogger(__name__)


class GeneInteractions(Gene):

    ''' GeneInteractions class define functions for building interations index type within gene index

    The interations index type is currently built by parsing the following:
    1. Refer section [INTACT] in download.ini for source files
    2. Refer section [BIOPLEX] in download.ini for source files

    Note: Most of the interaction data sources stores the interactions as binary interactions
    GeneA     GeneB
    100       728378
    100       345651
    645121    3312
    645121    55132
    645121    1020

    These files are parsed and entrezids are converted to ensemblids where needed.
    The interactors are grouped/clustered as below

    Grouping/clustering:
    100 => [728378, 345651]
    645121 => [3312, 55132, 1020]

    Final JSON structure that will be loaded
    {"interaction_source": "bioplex", "interactors": [{"interactor": "ENSG00000143416"},
                                                  {"interactor": "ENSG00000102043"},
                                                  {"interactor": "ENSG00000079387"},
                                                  {"interactor": "ENSG00000187231"}],
                                                  "_parent": "ENSG00000152213"}
    '''

    @classmethod
    def gene_interaction_parse(cls, download_file, stage_output_file, section, config=None):
        '''Function to delegate parsing of gene interaction files based on the file formats eg: psimitab'''
        if str(section._name) == "INTACT":
            cls._psimitab(download_file, stage_output_file, section, config)
        if str(section._name) == "BIOPLEX":
            cls._process_bioplex(download_file, stage_output_file, section, config)

    @classmethod
    def _process_bioplex(cls, download_file, stage_output_file, section, config):
        '''Function to process bioplex data files. Interactors are in first two columns, they are converted to
        ensembl ids and stored in temperory.out files
        Input File format:
        GeneA    GeneB    UniprotA    UniprotB    SymbolA    SymbolB    pW    pNI    pInt
        100    728378    P00813    A5A3E0    ADA    POTEF    2.38086E-09    0.000331856    0.999668142
        100    345651    P00813    Q562R1    ADA    ACTBL2    9.79E-18    0.211914437    0.788085563

        Output file format:
        interactorA    interactorB
        ENSG00000196839    ENSG00000196604
        ENSG00000196839    ENSG00000169067
        '''
        stage_output_file_handler = open(stage_output_file, 'w')
        mapped_counter = 0
        unmapped_ids = []
        stage_output_file_handler.write('interactorA' + '\t' + 'interactorB\n')

        gene_sets = []
        with open(download_file, encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader:
                gene_sets.extend([row['GeneA'], row['GeneB']])
        csvfile.close()

        ens_look_up = Gene._entrez_ensembl_lookup(gene_sets, section, config)

        with open(download_file, encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader:
                interactor_a = row['GeneA']
                interactor_b = row['GeneB']
                if interactor_a in ens_look_up and interactor_b in ens_look_up:
                    line = ens_look_up[interactor_a] + '\t' + ens_look_up[interactor_b] + '\n'
                    stage_output_file_handler.write(line)
                    mapped_counter += 1
                else:
                    line = interactor_a + '\t' + interactor_b + '\n'
                    unmapped_ids.append(interactor_a)
                    unmapped_ids.append(interactor_b)

        logger.debug("\n".join(unmapped_ids))
        logger.debug("Mapped {}  Unmapped {} " . format(mapped_counter, len(unmapped_ids)))

        stage_output_file_handler.close()
        cls._process_interaction_out_file(stage_output_file, section, False)

    @classmethod
    def _psimitab(cls, download_file, stage_output_file, section, config):
        '''Function to process intact psimitab data files
        Input file is the psimitab file

        Output file is:
        interactorA    interactorB    pubmed
        ENSG00000078053    ENSG00000159082    10542231
        ENSG00000078053    ENSG00000159082    10542231
        ENSG00000078053    ENSG00000159082    10542231
        '''
        abs_path_download_dir = os.path.dirname(download_file)
        zf = zipfile.ZipFile(download_file, 'r')

        import_file_exists = False

        if import_file_exists is False:
            stage_output_file_handler = open(stage_output_file, 'w')
            header_line = 'interactorA' + '\t' + 'interactorB' + '\t' + 'pubmed' + '\n'

            stage_output_file_handler.write(header_line)

            if 'intact.txt' in zf.namelist():

                print('Extracting the zip file...')
                target_path = zf.extract(member='intact.txt', path=abs_path_download_dir)
                line_number = 0
                with open(target_path, encoding='utf-8') as csvfile:
                    reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
                    for row in reader:
                            is_human_A = cls._check_tax_id(row['Taxid interactor A'], 'taxid:9606')
                            is_human_B = cls._check_tax_id(row['Taxid interactor B'], 'taxid:9606')

                            if is_human_A and is_human_B:
                                pass
                            else:
                                continue

                            # check for pubmed id/evidence
                            cleaned_pubmed_id = cls._clean_id(row['Publication Identifier(s)'], 'pubmed:\d+')

                            if cleaned_pubmed_id is None:
                                cleaned_pubmed_id = ''
                            # xref id
                            cleaned_xref_id_A = cls._clean_id(row['Xref(s) interactor A'], 'ensembl:ENSG\d+')
                            cleaned_xref_id_B = cls._clean_id(row['Xref(s) interactor B'], 'ensembl:ENSG\d+')

                            if (cleaned_xref_id_A is not None and
                               cleaned_xref_id_B is not None and
                               cleaned_xref_id_A != cleaned_xref_id_B):
                                line_number += 1
                                line = cleaned_xref_id_A + '\t' + cleaned_xref_id_B + '\t' + cleaned_pubmed_id + '\n'
                                stage_output_file_handler.write(line)

                stage_output_file_handler.close()
                cls._process_interaction_out_file(stage_output_file, section)
        else:
            cls._process_interaction_out_file(stage_output_file, section)

    @classmethod
    def _process_interaction_out_file(cls, target_path, section, include_evidence=True):
        '''Process the tab limited interaction output file to groups/cluster the interactors
        input file format:
        interactorA    interactorB    pubmed
        ENSG00000078053    ENSG00000159082    10542231
        ENSG00000078053    ENSG00000159082    10542231
        '''
        line_number = 0
        gene_interactors_dict = dict()
        evidence_id = None

        with open(target_path) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader:
                line_number += 1
                sys.stdout.write('.')
                # print(row)
                cleaned_xref_id_A = row['interactorA']
                cleaned_xref_id_B = row['interactorB']
                if include_evidence:
                    cleaned_pubmed_id = row['pubmed']
                    evidence_id = cleaned_pubmed_id

                if cleaned_xref_id_A == cleaned_xref_id_B:
                    continue

                interactorA = cleaned_xref_id_A
                interactorB = cleaned_xref_id_B

                # pass the interactors and get back the list
                (gene_interactors_list_a, gene_interactors_list_b) = cls._check_binary_interactions(gene_interactors_dict,  # @IgnorePep8
                                                                                                    interactorA,
                                                                                                    interactorB,
                                                                                                    evidence_id)
                gene_interactors_dict[interactorA] = gene_interactors_list_a
                gene_interactors_dict[interactorB] = gene_interactors_list_b

            cls._create_json_output_interaction(gene_interactors_dict, target_path, section)
            print('GENE INTERACTION STAGE COMPLETE')

    @classmethod
    def _create_json_output_interaction(cls, dict_container, target_file_path, section):
        '''Stores the output from _process_interaction_out_file function into JSON file'''
        count = 0
        dict_keys = dict_container.keys()
        json_target_file_path = target_file_path.replace(".out", ".json")
        interaction_source = section['source'].lower()

        load_mapping = True
        with open(json_target_file_path, mode='w', encoding='utf-8') as f:
            f.write('{"docs":[\n')

            for gene in dict_container:
                # decorate the list
                gene_list = dict_container[gene]
                list2json = cls.interaction_json_decorator(interaction_source, gene, gene_list)
                f.write(list2json)
                count += 1

                if len(dict_keys) == count:
                    f.write('\n')
                else:
                    f.write(',\n')

            f.write('\n]}')
        logger.debug("No. genes to load "+str(count))
        logger.debug("Json written to " + json_target_file_path)
        logger.debug("Load mappings")

        if load_mapping:
            status = cls._load_interaction_mappings(section)
            logger.debug(str(status))

    @classmethod
    def _load_interaction_mappings(cls, section):
        '''Load the mappings for interactions index type'''
        interaction_mapping = MappingProperties("interactions", "gene")
        interaction_mapping.add_property("interactors", "object")
        interaction_mapping.add_property("interaction_source", "string")
        load = Loader()
        idx = section['index']
        options = {"indexName": idx, "shards": 1}
        status = load.mapping(interaction_mapping, "interactions", analyzer=Loader.KEYWORD_ANALYZER, **options)
        return status

    @classmethod
    def _check_tax_id(cls, search_str, idpattern):
        '''Utility function to checks if the Taxid column matches the given taxid string'''
        p = re.compile(idpattern)
        m = p.search(search_str)

        if m:
            return True
        else:
            return False

    @classmethod
    def _clean_id(cls, search_str, idpattern):
        '''Utility function to split the given string and get the value alone
        eg: pubmed:\d+ returns the pubmed id alone'''
        p = re.compile(idpattern)
        m = p.search(search_str)
        matched_id = ''
        if m:
            matched_id = m.group()
            split_id = matched_id.split(sep=':')
            if split_id:
                return split_id[1]
        else:
            pass

        return None

    @classmethod
    def _group_binary_interactions(cls, binary_interactions=None):
        '''Function to group and expand binary interactions...
        Takes a list of binary interactions as argument and delegates to _check_binary_interactions for each pair'''
        gene_interactors_dict = dict()

        for i, j in binary_interactions:
            interactorA = str(i)
            interactorB = str(j)

            if interactorA == interactorB:
                continue

            (gene_interactors_list_i, gene_interactors_list_j) = cls._check_binary_interactions(gene_interactors_dict,
                                                                                                interactorA,
                                                                                                interactorB)
            gene_interactors_dict[i] = gene_interactors_list_i
            gene_interactors_dict[j] = gene_interactors_list_j

        return gene_interactors_dict

    @classmethod
    def _check_binary_interactions(cls, gene_interactors, interactorA, interactorB, evidence_id=None):
        '''Function to check if the interactors exists in the given dict...
        if present append to the existing list and if not add them as new list'''
        interactorA = str(interactorA)
        interactorB = str(interactorB)

        if interactorA == interactorB:
            return (gene_interactors[interactorA], gene_interactors[interactorB])

        existing_listA = None
        if interactorA in gene_interactors:
            existing_listA = gene_interactors[interactorA]

            if not any(interactorB in d for d in existing_listA):
                interactorB_evidence = {interactorB: evidence_id}
                existing_listA.append(interactorB_evidence)
                gene_interactors[interactorA] = existing_listA
        else:
            gene_interactors[interactorA] = [{interactorB: evidence_id}]

        existing_listB = None
        if interactorB in gene_interactors:
            existing_listB = gene_interactors[interactorB]

            if not any(interactorA in d for d in existing_listB):
                interactorA_evidence = {interactorA: evidence_id}
                existing_listB.append(interactorA_evidence)
                gene_interactors[interactorB] = existing_listB

        else:
            gene_interactors[interactorB] = [{interactorA: evidence_id}]

        return (gene_interactors[interactorA], gene_interactors[interactorB])

    @classmethod
    def interactor_json_decorator(cls, gene_interactor, evidence_key="pubmed"):
        '''Given a dict  {geneA:12345}, returns back formatted json string as
        {'interactor': 'geneA'} or {'interactor': 'geneA', 'pubmed':12345}
        '''
        interactor, evidence_value = list(gene_interactor.items())[0]

        if evidence_value is not None:
            evidence_value = str(evidence_value)
            json_str = {"interactor": interactor, evidence_key: evidence_value}
        else:
            interactor = str(interactor)
            json_str = {"interactor": interactor}

        return json_str

    @classmethod
    def interaction_json_decorator(cls, interaction_source, parent, gene_list):
        '''
        Given a interactor list, returns the json formatted interaction
         {"interaction_source": "bioplex", "interactors": [{"interactor": "ENSG00000143416"},
                                                  {"interactor": "ENSG00000102043"},
                                                  {"interactor": "ENSG00000079387"},
                                                  {"interactor": "ENSG00000187231"}],
                                                  "_parent": "ENSG00000152213"}
        {"interaction_source": "bioplex", "interactors": [{"interactor": "ENSG00000143416", "pubmed":"1234"},
                                                  {"interactor": "ENSG00000102043", "pubmed":"1234"},
                                                  {"interactor": "ENSG00000079387", "pubmed":"3456"},
                                                  {"interactor": "ENSG00000187231", "pubmed":"1234"}],
                                                  "_parent": "ENSG00000152213"}
        '''
        interactors = []
        parent = str(parent)
        for interactor in gene_list:
            interator_json = cls.interactor_json_decorator(interactor)
            interactors.append(interator_json)

        interaction_json_str = json.dumps({"interaction_source": interaction_source, "_parent": parent,
                                           "interactors": interactors}, sort_keys=True)

        return interaction_json_str
