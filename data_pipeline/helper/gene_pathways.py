import logging
from builtins import classmethod
import csv
import os
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
import json
from data_pipeline.helper.gene import Gene

logger = logging.getLogger(__name__)


class GenePathways(Gene):

    ''' GenePathways class defines functions for building pathway_genesets index type within gene index

    The pathway_genesets index type is currently built by parsing the following:
    1. Refer section [MSIGDB] in download.ini for source files
    '''
    @classmethod
    def gene_pathway_parse(cls, download_files, stage_output_file, section):
        ''' Function to delegate parsing of gene pathway files based on the file formats eg: gmt - genematrix  '''
        cls._genematrix(download_files, stage_output_file, section)

    @classmethod
    def _genematrix(cls, download_files, stage_output_file, section):
        '''Function to delegate parsing of pathway files based on the source eg: kegg, reactome, go'''
        abs_path_staging_dir = os.path.dirname(stage_output_file)
        source = None
        is_public = True
        for file in download_files:
            stage_output_file = abs_path_staging_dir + '/' + os.path.basename(file) + '.json'
            if 'kegg' in file:
                source = 'kegg'
            elif 'reactome' in file:
                source = 'reactome'
            elif 'biocarta' in file:
                source = 'biocarta'
            elif 'all' in file:
                source = 'GO'
            else:
                source = 'unknown'
            cls._process_pathway(file, stage_output_file, section, source, is_public)

    @classmethod
    def _process_pathway(cls, download_file, stage_output_file, section, source, is_public):
        '''Function to parse the pathway input files eg: kegg, reactome, go
        INPUT file format:
        Pathway name \t Pathyway url \t List of entrez ids
        REACTOME_RNA_POL_I_TRANSCRIPTION_TERMINATION
        http://www.broadinstitute.org/gsea/msigdb/cards/REACTOME_RNA_POL_I_TRANSCRIPTION_TERMINATION1022
        2068    2071    25885    284119    2965    2966    2967    2968    4331

        The entrez ids are converted to ensembl ids and logs are written to track the conversion rates (LESS/MORE/EQUAL)

        '''
        pathway_name = None
        pathway_url = None
        gene_sets = None

        convert2ensembl = True

        json_target_file_path = stage_output_file.replace(".out", ".json")
        log_target_file_path = json_target_file_path + ".log"

        json_target_file = open(json_target_file_path, mode='w', encoding='utf-8')
        json_target_file.write('{"docs":[\n')

        log_target_file = open(log_target_file_path, mode='w', encoding='utf-8')

        count = 0
        tmp_row_count_file = open(download_file, encoding='utf-8')
        row_count = sum(1 for row in tmp_row_count_file)
        logger.debug('Number of lines in the file ' + str(row_count))

        load_mapping = True
        with open(download_file, encoding='utf-8') as csvfile:
                    reader = csv.reader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)

                    for row in reader:
                        path_object = dict()
                        pathway_name = row[0]
                        pathway_url = row[1]
                        gene_sets = row[2:]
                        gene_sets_count = len(gene_sets)

                        converted_genesets = 0
                        if convert2ensembl:
                            converted_genesets = super()._convert_entrezid2ensembl(gene_sets, section)
                            converted_genesets_count = len(converted_genesets)
                            gene_sets = converted_genesets

                        less_more = gene_sets_count/converted_genesets_count
                        if less_more > 1:
                            diff = gene_sets_count - converted_genesets_count
                            diff_text = "Less("+str(diff) + ")"
                        elif less_more < 1:
                            diff = converted_genesets_count - gene_sets_count
                            diff_text = "More("+str(diff) + ")"
                        else:
                            diff = converted_genesets_count - gene_sets_count
                            diff_text = "Equal("+str(diff) + ")"

                        logger.debug(pathway_name + "\t" + str(gene_sets_count) + '/' + str(converted_genesets_count)
                                     + "\t" + diff_text + "\n")
                        log_target_file.write(pathway_name + "\t" + str(gene_sets_count) + '/' +
                                              str(converted_genesets_count) + "\t" + diff_text + "\n")

                        path_object["pathway_name"] = pathway_name
                        path_object["pathway_url"] = pathway_url
                        path_object["gene_sets"] = gene_sets
                        path_object["source"] = source
                        path_object["is_public"] = is_public
                        json_target_file.write(json.dumps(path_object))
                        count += 1
                        if row_count == count:
                            json_target_file.write('\n')
                        else:
                            json_target_file.write(',\n')

                    json_target_file.write('\n]}')

        logger.debug("No. genes to load "+str(count))
        logger.debug("Json written to " + json_target_file_path)
        logger.debug("Load mappings")

        if load_mapping:
            status = cls._load_pathway_mappings(section)
            print(status)

    @classmethod
    def _load_pathway_mappings(cls, section):
        '''Function to load the elastic mappings'''
        idx = section['index']
        idx_type = section['index_type']
        pathway_mapping = MappingProperties(idx_type)
        pathway_mapping.add_property("pathway_name", "string")
        pathway_mapping.add_property("pathway_url", "string")
        pathway_mapping.add_property("gene_sets", "string")
        pathway_mapping.add_property("source", "string")
        pathway_mapping.add_property("is_public", "string")
        load = Loader()
        options = {"indexName": idx, "shards": 1}
        status = load.mapping(pathway_mapping, idx_type, **options)
        return status