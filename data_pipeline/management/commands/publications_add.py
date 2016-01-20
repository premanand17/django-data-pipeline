''' Add publications from a list of PMIDs in a file. '''
from django.core.management.base import BaseCommand
from django.core.management import call_command
from data_pipeline import utils
from data_pipeline.helper.pubs import Pubs


class Command(BaseCommand):
    help = "Add a list of PMIDs from a file."

    def add_arguments(self, parser):
        parser.add_argument('--pmids',
                            dest='pmids',
                            help='File of PMIDs', required=True)
        parser.add_argument('--idx',
                            dest='idx',
                            help='Publications Index', required=True)

    def handle(self, *args, **options):
        file_of_pmids = options['pmids']
        idx = options['idx']
        pmids = [line.replace("'", "").strip() for line in open(file_of_pmids) if line.strip() != ""]

        print("Load PMIDs from file "+file_of_pmids)
        new_pmids = utils.PostProcess.get_new_pmids(pmids, idx)
        print("Number of new PMIDs = "+str(len(new_pmids)))

        stage_file = file_of_pmids + ".xml"
        print('Staging file: '+stage_file)
        Pubs.fetch_details(new_pmids, stage_file)
        call_command('index_search', indexType='publication', indexJson=stage_file, indexName=idx)
