''' Command line tool to get download sources and versions. '''
from django.core.management.base import BaseCommand, CommandError
import logging
from data_pipeline.utils import IniParser
import json

# Get an instance of a logger
logger = logging.getLogger(__name__)


class Ini2Json(IniParser):

    res = {}

    def convert(self, ini_file):
        ''' Covert the ini file. '''
        self.process_sections(self.read_ini(ini_file), ".")
        print(json.dumps(Ini2Json.res))

    def process_section(self, fname, section_dir_name, base_dir_path,
                        dir_path='.', section=None, stage='download', config=None):
        ''' Overrides L{IniParser.process_section} to process a section
        in the config file '''
        s = Ini2Json.res

        d = ''
        if 'location' in section:
            d += section['location']
            if 'files' in section:
                d += section['files']

        if d != '':
            s[section.name] = {'download': d}


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--ini',
                            dest='ini',
                            help='Input file defining downloads.')

    def handle(self, *args, **options):
        logger.debug(options)
        ini = Ini2Json()
        if options['ini']:
            if options['ini']:
                ini.convert(options['ini'])
            else:
                raise CommandError('--ini parameter not provided')
