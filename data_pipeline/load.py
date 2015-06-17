''' Loaders '''
import os
import logging
from .utils import IniParser
from django.core.management import call_command

# Get an instance of a logger
logger = logging.getLogger(__name__)


class IndexLoad(IniParser):
    ''' Load into Elastic '''

    def load(self, ini_file, dir_path, sections=None):
        ''' Download data defined in the ini file. '''
        return self.process_sections(self.read_ini(ini_file), dir_path, sections)

    def process_section(self, section_name, section_dir_name, base_dir_path, dir_path='.', section=None):
        ''' Overrides L{IniParser.process_section} to process a section
        in the config file '''
        if 'output' not in section:
            return
        stage_file = os.path.join(base_dir_path, 'STAGE', section_dir_name,
                                  section['output'] + '.json')

        logger.debug('Loading: '+stage_file + ' into ' + section['index'])
        call_command('index_search', indexType='auto', indexJson=stage_file, indexName=section['index'])
