from psm_utils.io import read_file
import traceback
import logging
import logging.config
import json
from datetime import datetime
from pyteomics import mgf, mzml
from urllib.parse import urlparse
from os.path import exists

RAW_FILE = r"blabla.raw"
IDENT_FILE = r"D:\dev\EuBIC-MS Developers meeting 2023\fragannot\data\Eubic2023\Human-Protein-Training_Trypsin.mzid"

class Parser:
    
    def __init__(self):
        # Set up logging
        logging.config.dictConfig(self.__load_log_config())
        self.logger = logging.getLogger(__name__)
        
    def read(self, raw_file, ident_file):
        try:
            #raw = 
            psm_list = read_file(ident_file)
        except Exception as ex:
            parser.logger.info(f"Couldn't read file. Exception:\n{traceback.format_exc()}")
        
    def __load_log_config(self):
        config = {}
        with open("log_conf.json", "r", encoding="utf-8") as fd:
            config = json.load(fd)
        return config
    
    # Check if the given filename is a url or not
    def __is_local(url):
        url_parsed = urlparse(url)
        if url_parsed.scheme in ('file', ''):
            return exists(url_parsed.path)
        return False

if __name__ == "__main__":
    parser = Parser()
    parser.read(RAW_FILE, IDENT_FILE)
    