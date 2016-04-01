''' Command line tool to get download sources and versions. '''
import datetime
import logging
from django.core.management.base import BaseCommand, CommandError
from data_pipeline.utils import IniParser


# Get an instance of a logger
logger = logging.getLogger(__name__)


class Ini2Template(IniParser):

    def convert(self, ini_file):
        ''' Covert the ini file. '''
        self.process_sections(self.read_ini(ini_file), ".")
        # print(json.dumps(Ini2Json.res))

    def process_section(self, fname, section_dir_name, base_dir_path,
                        dir_path='.', section=None, stage='download', config=None):
        ''' Overrides L{IniParser.process_section} to process a section
        in the config file '''
        d = ''

        if 'location' in section:
            d += section['location']
            if 'files' in section:
                d += section['files']

        if d != '':
            name = section['name'] if 'name' in section else section.name
            print('<div class="row">')
            print('<div class="col-md-2">' + name +
                  '</div><div class="col-md-8">' + d, end="</div>")

            if 'version' in section:
                print('<div class="col-md-2">'+section['version'], end="</div>")
            else:
                print('<div class="col-md-2">'+datetime.datetime.now().strftime("%d-%m-%y"), end="</div>")
            print("</div>")


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--ini',
                            dest='ini',
                            help='Input file defining downloads.')

    def handle(self, *args, **options):
        logger.debug(options)
        ini = Ini2Template()
        if options['ini']:
            if options['ini']:
                ini.convert(options['ini'])
            else:
                raise CommandError('--ini parameter not provided')
