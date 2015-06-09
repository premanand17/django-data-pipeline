''' Command line tool to manage downloads. '''
from django.core.management.base import BaseCommand, CommandError
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
        parser.add_argument('--ini',
                            dest='ini',
                            help='Input file defining downloads.')
        parser.add_argument('--download',
                            dest='download',
                            help='Download step',
                            action="store_true")

    def handle(self, *args, **options):

        if options['download']:
            if options['ini']:
                if not options['dir']:
                    raise CommandError('--dir parameter not provided')
                if Download().download_ini(options['ini'], options['dir']):
                    self.stdout.write("DOWNLOAD COMPLETE")
            else:
                if Download().download(options['url'], options['dir']):
                    self.stdout.write("DOWNLOAD COMPLETE")
