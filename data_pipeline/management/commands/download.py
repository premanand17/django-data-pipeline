''' Command line tool to manage downloads. '''
from django.core.management.base import BaseCommand
from data_pipeline.download import Download


class Command(BaseCommand):
    help = "Download data file(s)."

    def add_arguments(self, parser):
        parser.add_argument('url',
                            type=str,
                            help='Data URL')
        parser.add_argument('--dir',
                            dest='dir',
                            metavar="/download_path/",
                            help='Directory to store downloads.')

    def handle(self, *args, **options):

        Download().download(options['url'], options['dir'])
