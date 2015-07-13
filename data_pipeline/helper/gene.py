import logging
from builtins import classmethod
import json
from elastic.management.loaders.mapping import MappingProperties
from elastic.management.loaders.loader import Loader
import sys
import csv
import re
import zipfile
import os
logger = logging.getLogger(__name__)


class Gene(object):

    @classmethod
    def gene_info_parse(cls, gene_infos):
        gene_list = {}
        for gene_info in gene_infos:
            if gene_info.startswith('9606\t'):
                parts = gene_info.split('\t')
                gi = {}
                gi['symbol'] = parts[2]
                gi['synonyms'] = parts[4].split('|')
                gi['description'] = parts[8]
                gi['type'] = parts[9]
                gi['nomenclature_authority_symbol'] = parts[10]
                gi['full_name'] = parts[11]
                if parts[5] != '-':
                    cls._set_dbxrefs(parts[5], gi)
                gene_list[parts[1]] = gi
        return gene_list

    @classmethod
    def _set_dbxrefs(cls, dbxrefs, gi):
        arr = {}
        for dbxref in dbxrefs.split('|'):
            dbx = dbxref.rsplit(":", 1)
            if len(dbx) != 2:
                logger.warn('DBXREF PARSE: '+dbxref)
                continue
            if dbx[0].lower() == 'ensembl':
                if 'ensembl' in gi:
                    if not isinstance(gi['ensembl'], list):
                        gi['ensembl'] = [gi['ensembl']]
                        gi['ensembl'].append(dbx[1])
                    else:
                        gi['ensembl'] = dbx[1]
                else:
                    arr[dbx[0]] = dbx[1]
        gi['dbxrefs'] = arr

    @classmethod
    def gene_interaction_parse(cls, download_file, stage_output_file, section):
        cls._psimitab(download_file, stage_output_file, section)

    @classmethod
    def _psimitab(cls, download_file, stage_output_file, section):
        abs_path_download_dir = os.path.dirname(download_file)
        zf = zipfile.ZipFile(download_file, 'r')
        # print(zf.namelist())

        import_file_exists = False
        # target_path = '/dunwich/scratch/prem/tmp/download/intact.txt'

        if import_file_exists is not True:
            stage_output_file_handler = open(stage_output_file, 'w')
            header_line = 'interactorA' + '\t' + 'interactorB' + '\t' + 'pubmed' + '\n'

            stage_output_file_handler.write(header_line)

            if 'intact.txt' in zf.namelist():
                # base_dir_path = args[3]
                print('Extracting the zip file...')
                target_path = zf.extract(member='intact.txt', path=abs_path_download_dir)
                line_number = 0
                with open(target_path) as csvfile:
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
    def _process_interaction_out_file(cls, target_path, section):
        print('_process_interaction_out_file file called...start processing')

        dict_container = dict()
        line_number = 0

        with open(target_path) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            for row in reader:
                line_number += 1
                sys.stdout.write('.')
                # print(row)
                cleaned_xref_id_A = row['interactorA']
                cleaned_xref_id_B = row['interactorB']
                cleaned_pubmed_id = row['pubmed']

                if cleaned_xref_id_A == cleaned_xref_id_B:
                    continue

                dict_interactorA = {'interactor': cleaned_xref_id_A, 'pubmed': cleaned_pubmed_id}
                dict_interactorB = {'interactor': cleaned_xref_id_B, 'pubmed': cleaned_pubmed_id}

                if cleaned_xref_id_A in dict_container:
                    current_A = dict_container[cleaned_xref_id_A]
                    if dict_interactorB not in current_A:
                        current_A.append(dict_interactorB)
                        dict_container[cleaned_xref_id_A] = current_A
                else:
                    dict_container[cleaned_xref_id_A] = [dict_interactorB]

                if cleaned_xref_id_B in dict_container:
                    current_B = dict_container[cleaned_xref_id_B]
                    if dict_interactorA not in current_B:
                        current_B.append(dict_interactorA)
                        dict_container[cleaned_xref_id_B] = current_B
                else:
                    dict_container[cleaned_xref_id_B] = [dict_interactorA]

            cls._create_json_output(dict_container, target_path, section)
            print('GENE INTERACTION STAGE COMPLETE')

    @classmethod
    def _create_json_output(cls, dict_container, target_file_path, section):
        count = 0
        dict_keys = dict_container.keys()
        print('No of docs to load ' + str(len(dict_keys)))
        json_target_file_path = target_file_path.replace(".out", ".json")
        print(json_target_file_path)
        # json_target_file_path = target_file_path + '.json'

        load_mapping = True
        with open(json_target_file_path, mode='w', encoding='utf-8') as f:
            f.write('{"docs":[\n')

            for gene in dict_container:
                int_object = dict()
                int_object["_id"] = gene
                int_object["_parent"] = gene
                int_object["interactors"] = dict_container[gene]
                int_object["interaction_source"] = section['source']
                f.write(json.dumps(int_object))
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
            print(str(status))

    @classmethod
    def _load_interaction_mappings(cls, section):
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
        p = re.compile(idpattern)
        m = p.search(search_str)

        if m:
            return True
        else:
            return False

    @classmethod
    def _clean_id(cls, search_str, idpattern):
        p = re.compile(idpattern)
        m = p.search(search_str)
        matched_id = ''
        if m:
            # print('Match found: ', m.group())
            matched_id = m.group()
            split_id = matched_id.split(sep=':')
            if split_id:
                return split_id[1]
        else:
            pass

        return None
