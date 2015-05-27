''' Command line tool to manage downloads. '''
from django.core.management.base import BaseCommand
from data_pipeline.download import Download


class Command(BaseCommand):
    help = "Download data file(s)."

    def add_arguments(self, parser):
        parser.add_argument('--url',
                            type=str,
                            help='Data URL')
        parser.add_argument('--dir',
                            dest='dir',
                            metavar="/download_path/",
                            help='Directory to store downloads.')
        parser.add_argument('--input',
                            dest='input',
                            help='Input file defining downloads.')
        parser.add_argument('--ini',
                            dest='ini',
                            help='Input file defining downloads.')

    def handle(self, *args, **options):
        if options['input']:
            Download().download_file(options['input'], options['dir'])
        elif options['ini']:
            Download().download_ini(options['ini'], options['dir'])
        else:
            Download().download(options['url'], options['dir'])
