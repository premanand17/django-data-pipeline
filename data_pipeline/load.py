''' Loaders '''
import os
import logging
from .utils import IniParser
from django.core.management import call_command
from data_pipeline.utils import pre_process

# Get an instance of a logger
logger = logging.getLogger(__name__)


class IndexLoad(IniParser):
    ''' Load into Elastic '''

    def load(self, ini_file, dir_path, sections=None):
        ''' Download data defined in the ini file. '''
        return self.process_sections(self.read_ini(ini_file), dir_path, sections)

    @pre_process
    def process_section(self, section_name, section_dir_name, base_dir_path,
                        dir_path='.', section=None, stage='load'):
        ''' Overrides L{IniParser.process_section} to process a section
        in the config file '''
        if 'output' in section:
            stage_file = os.path.join(base_dir_path, 'STAGE', section_dir_name,
                                      section['output'] + '.json')
        elif 'files' in section:
            stage_file = os.path.join(base_dir_path, 'STAGE', section_dir_name,
                                      section['files'] + '.json')
        else:
            return

        if not os.path.exists(stage_file):
            logger.error('File does not exist: '+stage_file)
            return
        logger.debug('Loading: '+stage_file + ' into ' + section['index'])
        call_command('index_search', indexType=section['index_type'], indexJson=stage_file, indexName=section['index'])
