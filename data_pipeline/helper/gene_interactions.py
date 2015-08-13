import logging
from builtins import classmethod
import sys
import csv
import re
import zipfile
import os
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
import json
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
    def gene_interaction_parse(cls, download_file, stage_output_file, section):
        '''Function to delegate parsing of gene interaction files based on the file formats eg: psimitab'''
        if str(section._name) == "INTACT":
            cls._psimitab(download_file, stage_output_file, section)
        if str(section._name) == "BIOPLEX":
            cls._process_bioplex(download_file, stage_output_file, section)

    @classmethod
    def _process_bioplex(cls, download_file, stage_output_file, section):
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
        unmapped_counter = 0
        unmapped_ids = []
        header_line = 'interactorA' + '\t' + 'interactorB\n'
        stage_output_file_handler.write(header_line)

        log_target_file = stage_output_file + ".log"
        log_target_file_handler = open(log_target_file, mode='w', encoding='utf-8')
        with open(download_file, encoding='utf-8') as csvfile:
                    reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
                    for row in reader:
                        gene_sets = []
                        interactor_a = row['GeneA']
                        interactor_b = row['GeneB']
                        gene_sets.append(interactor_a)
                        gene_sets.append(interactor_b)

                        ensembl_ids = super()._convert_entrezid2ensembl(gene_sets, section,
                                                                        log_target_file_handler, True)
                        if(len(ensembl_ids) == 2):
                            line = ensembl_ids[0] + '\t' + ensembl_ids[1] + '\n'
                            stage_output_file_handler.write(line)
                            mapped_counter += 1
                        else:
                            line = interactor_a + '\t' + interactor_b + '\n'
                            unmapped_counter += 1
                            unmapped_ids.append(interactor_a)
                            unmapped_ids.append(interactor_b)

        logger.debug("\n".join(unmapped_ids))
        logger.debug("Mapped {}  Unmapped {} " . format(mapped_counter, unmapped_counter))

        stage_output_file_handler.close()
        # log_target_file_handler.close()
        cls._process_interaction_out_file(stage_output_file, section, False)

    @classmethod
    def _psimitab(cls, download_file, stage_output_file, section):
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

        if import_file_exists is not True:
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
                            # print(row)
                            # check for taxid
                            re.compile('taxid:9606')
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
        dict_container = dict()
        line_number = 0
        interaction_source = section['source'].lower()
        gene_interactors_dict = dict()
        evidence_dict = dict()

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
                if include_evidence:
                    (gene_interactors_list_a, gene_interactors_list_b, evidence_list_a, evidence_list_b) = cls._check_binary_interactions(gene_interactors_dict,  # @IgnorePep8
                                                                                                                                          interactorA, # @IgnorePep8
                                                                                                                                          interactorB, # @IgnorePep8
                                                                                                                                          evidence_dict, # @IgnorePep8
                                                                                                                                          evidence_id) # @IgnorePep8
                else:
                    (gene_interactors_list_a, gene_interactors_list_b, evidence_list_a, evidence_list_b) = cls._check_binary_interactions(gene_interactors_dict, # @IgnorePep8
                                                                                                                                          interactorA, # @IgnorePep8
                                                                                                                                          interactorB) # @IgnorePep8

                gene_interactors_dict[interactorA] = gene_interactors_list_a
                gene_interactors_dict[interactorB] = gene_interactors_list_b

                if include_evidence:
                    evidence_key = "pubmed"
                    json_interaction_a = cls.interaction_json_decorator(interaction_source, interactorA,
                                                                        gene_interactors_list_a,
                                                                        evidence_key, evidence_list_a)
                    json_interaction_b = cls.interaction_json_decorator(interaction_source,
                                                                        interactorB,
                                                                        gene_interactors_list_b,
                                                                        evidence_key,
                                                                        evidence_list_b)
                else:
                    json_interaction_a = cls.interaction_json_decorator(interaction_source, interactorA,
                                                                        gene_interactors_list_a)
                    json_interaction_b = cls.interaction_json_decorator(interaction_source, interactorB,
                                                                        gene_interactors_list_b)

                dict_container[interactorA] = json_interaction_a
                dict_container[interactorB] = json_interaction_b

            cls._create_json_output_interaction(dict_container, target_path, section)
            print('GENE INTERACTION STAGE COMPLETE')

    @classmethod
    def _create_json_output_interaction(cls, dict_container, target_file_path, section):
        '''Stores the output from _process_interaction_out_file function into JSON file'''
        count = 0
        dict_keys = dict_container.keys()
        json_target_file_path = target_file_path.replace(".out", ".json")

        load_mapping = False
        with open(json_target_file_path, mode='w', encoding='utf-8') as f:
            f.write('{"docs":[\n')

            for gene in dict_container:
                f.write(dict_container[gene])
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
        interaction_mapping.add_property("interactors", "object", index="not_analyzed")
        interaction_mapping.add_property("interaction_source", "string")
        load = Loader()
        idx = section['index']
        options = {"indexName": idx, "shards": 1}
        status = load.mapping(interaction_mapping, "interactions", **options)
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
            i = str(i)
            j = str(j)

            if i == j:
                continue

            (gene_interactors_list_i, gene_interactors_list_j, evidence_list_a, evidence_list_b) = cls._check_binary_interactions(gene_interactors_dict,  # @IgnorePep8 @UnusedVariable
                                                                                                                                  i, j)   # @IgnorePep8
            gene_interactors_dict[i] = gene_interactors_list_i
            gene_interactors_dict[j] = gene_interactors_list_j

        return gene_interactors_dict

    @classmethod
    def _check_binary_interactions(cls, gene_interactors, interactorA, interactorB, evidence_dict={}, evidence_id=None):
        '''Function to check if the interactors exists in the given dict...
        if present append to the existing list and if not add them as new list'''
        i = str(interactorA)
        j = str(interactorB)

        if i == j:
            return (gene_interactors[i], gene_interactors[j])

        existing_list = None
        evidence_existing_list = None
        if i in gene_interactors:
            existing_list = gene_interactors[i]

            if j not in existing_list:
                existing_list.append(j)
                gene_interactors[i] = existing_list

            if i in evidence_dict:
                evidence_existing_list = evidence_dict[i]

                if evidence_existing_list and evidence_id is not None:
                    evidence_existing_list.append(evidence_id)
                    evidence_dict[i] = evidence_existing_list
        else:
            gene_interactors[i] = [j]
            if evidence_id is not None:
                evidence_dict[i] = [evidence_id]

        existing_list = None
        evidence_existing_list = None
        if j in gene_interactors:
            existing_list = gene_interactors[j]

            if i not in existing_list:
                existing_list.append(i)
                gene_interactors[j] = existing_list

            if j in evidence_dict:
                evidence_existing_list = evidence_dict[j]

                if evidence_existing_list and evidence_id is not None:
                    evidence_existing_list.append(evidence_id)
                    evidence_dict[j] = evidence_existing_list
        else:
            gene_interactors[j] = [i]
            if evidence_id is not None:
                evidence_dict[j] = [evidence_id]

        if evidence_dict and len(evidence_dict) > 0:
            return (gene_interactors[i], gene_interactors[j], evidence_dict[i], evidence_dict[j])
        else:
            return (gene_interactors[i], gene_interactors[j], [], [])

    @classmethod
    def interactor_json_decorator(cls, gene_interactor, evidence_key=None, evidence_value=None):

        if evidence_key:
            json_str = {"interactor": gene_interactor, evidence_key: evidence_value}
        else:
            json_str = {"interactor": gene_interactor}

        return json_str

    @classmethod
    def interaction_json_decorator(cls, interaction_source, parent, gene_list, evidence_key=None, evidence_list=None):  # @IgnorePep8
        '''
         {"interaction_source": "bioplex", "interactors": [{"interactor": "ENSG00000143416"},
                                                  {"interactor": "ENSG00000102043"},
                                                  {"interactor": "ENSG00000079387"},
                                                  {"interactor": "ENSG00000187231"}],
                                                  "_parent": "ENSG00000152213"}
        '''
        interactors = []
        if evidence_list:
            for interactor, evidence_value in zip(gene_list, evidence_list):
                interator_json = cls.interactor_json_decorator(interactor, evidence_key, evidence_value)
                interactors.append(interator_json)
        else:
            for interactor in gene_list:
                interator_json = cls.interactor_json_decorator(interactor)
                interactors.append(interator_json)

        interaction_json_str = json.dumps({"interaction_source": interaction_source, "_parent": parent,
                                           "interactors": interactors}, sort_keys=True)

        return interaction_json_str
