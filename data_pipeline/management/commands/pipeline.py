''' Command line tool to manage downloads. '''
from django.core.management.base import BaseCommand, CommandError
from data_pipeline.download import Download
from data_pipeline.stage import Stage
from data_pipeline.load import IndexLoad
import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    ''' Command lines for downloading and loading data. The order for the gene index below is important
    for the correct mappings to be loaded. Elasticsearch (2.0.0) requires the mapping for interactions
    to be in place before the parent (gene) objects - so INTACT is loaded before the main gene types.

    Gene History:
    ./manage.py pipeline --dir tmp --ini download.ini --sections GENE_HISTORY --steps download load

    Gene Interactions:
    ./manage.py pipeline --dir tmp --ini download.ini --sections INTACT --steps download stage load

    Gene:
    ./manage.py pipeline --dir tmp --ini download.ini --sections ENSEMBL_GENE --steps download stage load
    ./manage.py pipeline --dir tmp --ini download.ini --sections GENE2ENSEMBL --steps download load
    ./manage.py pipeline --dir tmp --ini download.ini --sections ENSMART_GENE --steps download load
    ./manage.py pipeline --dir tmp --ini download.ini --sections GENE_INFO --steps download load
    ./manage.py pipeline --dir tmp --ini download.ini --sections GENE_PUBS --steps download load
    ./manage.py pipeline --dir tmp --ini download.ini --sections ENSMART_HOMOLOG --steps download load
    ./manage.py pipeline --dir tmp --ini download.ini --sections ENSEMBL2MGI --steps download load

    Gene Interactions:
    ./manage.py pipeline --dir tmp --ini download.ini --sections BIOPLEX --steps download stage load

    Gene Pathways/Genesets:
    ./manage.py pipeline --dir tmp --ini download.ini --sections MSIGDB --steps download stage load

    Update gene suggester weighting:
    python criteria_suggester.py gene


    Marker:
    ./manage.py pipeline --dir /dbSNP/human/144/ --ini download.ini --sections DBSNP --steps download load
    ./manage.py pipeline --dir /dbSNP/human/144/ --ini download.ini --sections RSMERGEARCH --steps download load

    ./manage.py pipeline --dir tmp --ini download.ini  --sections IMMUNOCHIP_MYSQL --steps load

    '''
    help = "Download data file(s)"

    def add_arguments(self, parser):
        parser.add_argument('--url',
                            type=str,
                            help='Data URL')
        parser.add_argument('--dir',
                            dest='dir',
                            metavar="/download_path/",
                            help='Directory to store downloads.')
        parser.add_argument('--ini',
                            dest='ini',
                            help='Input file defining downloads.')
        parser.add_argument('--sections',
                            dest='sections',
                            help='Comma separated section names (e.g. BANDS) to download [default: all].')
        parser.add_argument('--steps',
                            dest='steps',
                            help='Steps to run [download load]',
                            nargs='+', required=True)

    def handle(self, *args, **options):
        logger.debug(options)
        if 'download' in options['steps']:
            if options['ini']:
                if not options['dir']:
                    raise CommandError('--dir parameter not provided')
                if Download().download_ini(options['ini'], options['dir'], options['sections']):
                    self.stdout.write("DOWNLOAD COMPLETE")
            else:
                if Download().download(options['url'], options['dir']):
                    self.stdout.write("DOWNLOAD COMPLETE")
        if 'stage' in options['steps']:
            Stage().stage(options['ini'], options['dir'], options['sections'])
        if 'load' in options['steps']:
            IndexLoad().load(options['ini'], options['dir'], options['sections'])
