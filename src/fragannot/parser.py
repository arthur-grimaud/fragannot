from psm_utils.io import read_file
from psm_utils import Peptidoform, PSM
import traceback
import logging
import logging.config
import json
import requests
from datetime import datetime
from pyteomics import mgf, mzml, mzid
from urllib.parse import urlparse
from os.path import exists, splitext, split
from urllib.request import urlopen
import os
import re
import gzip

RAW_FILE = r"../../data/2020_09_90_092320_Hazbun_SigmaGluC_CID_orbiorbi.mgf"
IDENT_FILE = r"../../data/2020_09_90_092320_Hazbun_Sigma_GluC_CID_orbiorbi_peptide.mzid"

#RAW_FILE = r"https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/08/PXD006552/generated/Pt1_F_100K_tech-rep1.pride.mgf.gz"
#IDENT_FILE = r"https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/08/PXD006552/peptides_1_1_0.mzid.gz"

class Parser:
    def __init__(self):
        # Set up logging
        logging.config.dictConfig(self.__load_log_config())
        self.logger = logging.getLogger(__name__)
        self.spectra = None
        self.psm_list = None

    def read(self, raw_file, ident_file, file_format):
        """ Read and process raw file and identification file.
        
        Parameters
        ----------
        raw_file : str
            Path or url to the raw file
        ident_file_path : str
            Path or url to the identification file
        """
        try:
            if self.__is_url(raw_file):
                self.logger.info("Raw file is not local, try to download.")
                raw_file_name = self.__get_file_from_url(raw_file)
            else:
                if os.path.exists(raw_file):
                    raw_file_name = raw_file
                else:
                    self.logger.error(f"File doesn't exist: {raw_file}")
                    raise Exception("File doesn't exist")
            if self.__is_url(ident_file):
                self.logger.info("Ident file is not local, try to download.")
                ident_file_name = self.__get_file_from_url(ident_file)
            else:
                if os.path.exists(ident_file):
                    ident_file_name = ident_file
                else:
                    self.logger.error(f"File doesn't exist: {raw_file}")
                    raise Exception("File doesn't exist")
            self.__load(raw_file_name, ident_file_name, file_format)
        except Exception as ex:
            self.logger.error(f"Couldn't read file. Exception:\n{traceback.format_exc()}")

        self.logger.info(f"Read {len(self.psm_list)} PSMs from identification file")
        
        count = 0
        for psm in self.psm_list:
            try:
                spectrum = self.spectra.get_by_id(psm["spectrum_id"])
                psm.spectrum = {"mz": spectrum["m/z array"],
                                "intensity": spectrum["intensity array"]}
                count += 1
                if count % 256 == 0:
                    print(f'\r{count} spectra processed')
                    
            except KeyError:
                self.logger.warning(f'SpectrumId - {psm["spectrum_id"]} not found')
        output_fpath  = os.path.splitext(raw_file_name)[0] + '.json'
        self.output_fname = os.path.basename(output_fpath)
        return self.psm_list
        
    def __load(self, raw_file_path, ident_file_path, file_format):
        """ Load raw file and identification file.
        
        Parameters
        ----------
        raw_file_path : str
            Path to the raw file
        ident_file_path : str
            Path to the identification file
        """
        self.spectra = self.__read_raw_file(raw_file_path)
        self.psm_list = self.__read_id_file(ident_file_path, file_format)
        
    def __read_raw_file(self, file_path):
        """ Read raw file 
        
        Parameters
        ----------
        file_path : str
            Path to the raw file
        """
    
        # Infer filetype from file path
        extension = splitext(file_path)[1]

        if extension.lower() == ".mzml":
            self.logger.info(f"Inferred mzML format from {file_path}")
            return mzml.MzML(file_path)
        elif extension.lower() == ".mgf":
            self.logger.info(f"Inferred MGF format from {file_path}")
            return mgf.IndexedMGF(file_path)
        
        else:
            self.logger.error(
                f"Cannot infer format from {file_path}, only mzML and MGF formats are supported"
            )
            raise Exception("Unsupported spectra file format")
    
    def __read_id_file(self, file_path, file_format):
        """ Read identification file more generously then psm_utils
        
        Parameters
        ----------
        file_path : str
            Path to the raw file
        file_format : str
            Identification file format
        """
        extension = splitext(file_path)[1]
        
        if extension.lower() == ".mzid" or extension.lower() == ".mzidentml":
            result = []
            for psm in mzid.MzIdentML(file_path):
                spectrumID = psm['spectrumID']
                sequence = psm['SpectrumIdentificationItem'][0]['PeptideEvidenceRef'][0]['PeptideSequence']
                modifications = psm['SpectrumIdentificationItem'][0]['PeptideEvidenceRef'][0].get('Modification', None)
                filename = split(psm['location'])[1]
                if not modifications is None:
                    aas = [''] + [aa for aa in sequence] + ['']
                    for mod in modifications:
                        loc = mod['location']
                        mass = mod['monoisotopicMassDelta']
                        if 'residues' in mod.keys():
                            res = mod['residues'][0]
                            if loc > 0 and not res == aas[loc]:
                                raise Exception(f'Mismatch {modifications} {sequence} {res} {aas[loc]} {loc}')
                        try:        
                            aas[loc] += f'[+{mass}]' if mass > 0 else f'[{mass}]'
                        except Exception:
                            print(psm)
                            break
                
                    sequence = ''.join(aas[1:-1])
                    
                    if aas[0] != '':
                        sequence = f'{aas[0]}-{sequence}'
                    
                    if aas[-1] != '':
                        sequence = f'{sequence}-{aas[-1]}'
    
                result.append(PSM(peptidoform=Peptidoform(sequence), run=filename, spectrum_id=spectrumID))
        
            return result
        
        else:
            return read_file(file_path, filetype=file_format)

    def __load_log_config(self):
        """ Load log configurations. """
        config = {}
        with open("log_conf.json", "r", encoding="utf-8") as fd:
            config = json.load(fd)
        return config
        
    def __uncompress(self, file_path, block_size=65536):
        new_file_path = file_path
        if file_path.endswith('.gz'):
            new_file_path = file_path[:-3]
            with gzip.open(file_path, 'rb') as s_file, open(new_file_path, 'wb') as d_file:
                while True:
                    block = s_file.read(block_size)
                    if not block:
                        break
                    else:
                        d_file.write(block)
        return new_file_path

    # Check if the given filename is a url or not
    def __is_url(self, url):
        is_url = True
        url_parsed = urlparse(url)
        if url_parsed.scheme in ("file", ""):
            is_url = False
        return is_url
        
    def __get_file_from_url(self, url):
        """ Download file from url and return the path to the file
        
        Parameters
        ----------
        url : str
            URL to the file to download
        """
        if not os.path.exists(r".\downloads"):
            os.makedirs(r".\downloads")
        
        r = requests.get(url, allow_redirects=True)
        if url.find('/'):
            file_path = ".\downloads\\" + url.rsplit('/', 1)[1]
        else:
            file_path = ".\downloads\\" + self.__getFilename_fromCd(r.headers.get('content-disposition'))
                    
        open(file_path, 'wb').write(r.content)
        file_path = self.__uncompress(file_path)
        
        return file_path
        
    def __getFilename_fromCd(self, cd):
        """ Gets filename from content disposition
        
        Parameters
        ----------
        cd : str
            Content disposition
        """
        if not cd:
            return None
        fname = re.findall('filename=(.+)', cd)
        print(fname)
        if len(fname) == 0:
            return None
        return fname[0]


if __name__ == "__main__":
    parser = Parser()
    parser.read(RAW_FILE, IDENT_FILE)
    print(parser.output_fname)
