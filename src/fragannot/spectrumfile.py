import re
from pyteomics import mzml, mgf
from os.path import splitext
from collections import defaultdict

class SpectrumFile:
    def __init__(self, file_path):
        self.indices = defaultdict(dict)
        self.spectra_source = None
        self.file_format = None
        self._build_index = None
        
        self._load(file_path)
        self._build_index()
        
    def _load(self, file_path):
        extension = splitext(file_path)[1]

        if extension.lower() == ".mzml":
            print(f"Inferred mzML format from {file_path}")
            self.spectra_source = mzml.MzML(file_path)
            self.file_format = 'mzml'
            self._build_index = self._index_MZML
            
        elif extension.lower() == ".mgf":
            print(f"Inferred MGF format from {file_path}")
            self.spectra_source = mgf.IndexedMGF(file_path)
            self.file_format = 'mgf'
            self._build_index = self._index_MGF
        
        else:
            print(
                f"Cannot infer format from {file_path}, only mzML and MGF formats are supported"
            )
            raise Exception("Unsupported spectra file format")
    
    def get_by_id(self, id_string):
        return self.spectra_source.get_by_id(id_string)
    
    def _index_MGF(self):
        for index in range(len(self.spectra_source)):
            params = self.spectra_source[index]['params']
            if 'title' in params.keys():
                self.indices['title'][params['title']] = index
                
            if 'scans' in params.keys():
                self.indices['scan'][params['scans']] = index
    
    def _index_MZML(self):
        for spectrum in self.spectra_source:
            spectrumID = spectrum['id']
            index = spectrum['index']
            
            scan_match = re.match(r'.+scan(?:\s+)?=(?:\s+)?(\d+)', spectrumID)
            if not scan_match is None: 
                self.indices['scan'][scan_match.group(1)] = index
            
            self.indices['id'][spectrumID] = index            

if __name__ == '__main__':
    import sys
    spectra = SpectrumFile(sys.argv[1])
        