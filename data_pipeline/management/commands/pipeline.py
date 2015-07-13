''' Command line tool to manage downloads. '''
from django.core.management.base import BaseCommand, CommandError
from data_pipeline.download import Download
from data_pipeline.stage import Stage
from data_pipeline.load import IndexLoad


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
        parser.add_argument('--sections',
                            dest='sections',
                            help='Comma separated section names (e.g. BANDS) to download [default: all].')
        parser.add_argument('--steps',
                            dest='steps',
                            help='Steps to run [download load]',
                            nargs='+', required=True)

    def handle(self, *args, **options):
        print(options)

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
