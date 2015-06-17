''' Command line tool to manage downloads. '''
from django.core.management.base import BaseCommand, CommandError
from data_pipeline.download import Download
from data_pipeline.load import IndexLoad


class Command(BaseCommand):
    help = "Download data file(s)."

    def add_arguments(self, parser):
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
#         parser.add_argument('--download',
#                             dest='download',
#                             help='Download step',
#                             action="store_true")

    def handle(self, *args, **options):

        if not options['dir']:
            raise CommandError('--dir parameter not provided')
        if Download().download_ini(options['ini'], options['dir'], options['sections']):
            self.stdout.write("DOWNLOAD COMPLETE")
        IndexLoad().load(options['ini'], options['dir'], options['sections'])
