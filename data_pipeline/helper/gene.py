import logging
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
