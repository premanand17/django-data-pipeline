''' Staging data ready for loading '''
import os
import logging
from .utils import IniParser
from .utils import post_process

# Get an instance of a logger
logger = logging.getLogger(__name__)


class Stage(IniParser):
    ''' Load into Elastic '''

    def stage(self, ini_file, dir_path, sections=None):
        ''' Download data defined in the ini file. '''
        return self.process_sections(self.read_ini(ini_file), dir_path, sections)

    @post_process
    def process_section(self, section_name, section_dir_name, base_dir_path,
                        dir_path='.', section=None, stage='stage'):
        ''' Overrides L{IniParser.process_section} to process a section
        in the config file '''
        if 'output' not in section:
            return False
        download_file = os.path.join(base_dir_path, 'DOWNLOAD', section_dir_name,
                                     section['output'])

        if not os.path.exists(download_file):
            logger.error('File does not exist: '+download_file)
            return False
        logger.debug('Process: '+download_file)
        return True
