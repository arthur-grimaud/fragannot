from psm_utils.io import read_file
from psm_utils import Peptidoform, PSM
import traceback
import logging
import logging.config
import json
from datetime import datetime
from pyteomics import mgf, mzml, mzid
from urllib.parse import urlparse
from os.path import exists, splitext, split

RAW_FILE = r"../../data/2020_09_90_092320_Hazbun_SigmaGluC_CID_orbiorbi.mgf"
IDENT_FILE = r"../../data/2020_09_90_092320_Hazbun_Sigma_GluC_CID_orbiorbi_peptide.mzid"


class Parser:
    def __init__(self):
        # Set up logging
        logging.config.dictConfig(self.__load_log_config())
        self.logger = logging.getLogger(__name__)
        self.spectra = None
        self.psm_list = None

    def load(self, raw_file, ident_file):
        self.spectra = self.read_raw_file(raw_file)
        self.psm_list = self.read_id_file(ident_file)
            
    def read(self, raw_file, ident_file):
        try:
            self.load(raw_file, ident_file)
        except Exception as ex:
            self.logger.error(f"Couldn't read file. Exception:\n{traceback.format_exc()}")

        self.logger.info(f"Read {len(self.psm_list)} from identification file")
        
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
            
            
        return self.psm_list    
            
    def read_raw_file(self, file_path):
        # infer filetype from file path
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
    
    def read_id_file(self, file_path):
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
            return read_file(file_path)

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
