from psm_utils.io import read_file
import traceback
import logging
import logging.config
import json
from datetime import datetime
from pyteomics import mgf, mzml
from urllib.parse import urlparse
from os.path import exists, splitext

RAW_FILE = r"../../data/Human-Protein-Training_Trypsin.mzML"
IDENT_FILE = r"../../data/Human-Protein-Training_Trypsin.mzid"


class Parser:
    def __init__(self):
        # Set up logging
        logging.config.dictConfig(self.__load_log_config())
        self.logger = logging.getLogger(__name__)
        self.spectra = None
        self.psm_list = None

    def load(self, raw_file, ident_file):
        self.spectra = self.read_raw_file(raw_file)
        self.psm_list = read_file(ident_file)

    def read(self, raw_file, ident_file):
        try:
            self.load(raw_file, ident_file)
        except Exception as ex:
            self.logger.error(f"Couldn't read file. Exception:\n{traceback.format_exc()}")

        self.logger.info(f"Read {len(self.psm_list)} from identification file")

        for psm in self.psm_list:
            spectrum = self.spectra.get_by_id(psm["spectrum_id"])
            psm.spectrum = {"mz": spectrum["m/z array"], "intensity": spectrum["intensity array"]}

        return self.psm_list

    def read_raw_file(self, file_path):
        # infer filetype from file path
        extension = splitext(file_path)[1]

        if extension.lower() == ".mzml":
            self.logger.info(f"Inferred mzML format from {file_path}")
            return mzml.MzML(file_path)
        elif extension.lower() == ".mgf":
            self.logger.info(f"Inferred MGF format from {file_path}")
            return mgf.MGF(file_path)
        else:
            self.logger.error(
                f"Cannot infer format from {file_path}, only mzML and MGF formats are supported"
            )
            raise Exception("Unsupported spectra file format")

    def __load_log_config(self):
        config = {}
        with open("log_conf.json", "r", encoding="utf-8") as fd:
            config = json.load(fd)
        return config

    # Check if the given filename is a url or not
    def __is_local(url):
        url_parsed = urlparse(url)
        if url_parsed.scheme in ("file", ""):
            return exists(url_parsed.path)
        return False


if __name__ == "__main__":
    parser = Parser()
    parser.read(RAW_FILE, IDENT_FILE)
