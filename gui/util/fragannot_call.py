#!/usr/bin/env python3

from util.fragannot_codon import FragannotCodon

import os
import random
from datetime import datetime
from typing import Dict
from typing import List
from typing import BinaryIO

def fragannot_call(spectrum_file: BinaryIO,
                   identifications_file: BinaryIO,
                   tolerance: float,
                   fragment_types: List[str],
                   charges: List[str],
                   losses: List[str],
                   deisotope: bool,
                   file_format: str = "infer") -> Dict:

    output_name_prefix = datetime.now().strftime("%b-%d-%Y_%H-%M-%S") + "_" + str(random.randint(10000, 99999))

    # write uploaded files to tmp directory
    with open(output_name_prefix + spectrum_file.name, "wb") as f1:
        f1.write(spectrum_file.getbuffer())
    with open(output_name_prefix + identifications_file.name, "wb") as f2:
        f2.write(identifications_file.getbuffer())

    # run fragannot
    frag = FragannotCodon()
    fragannot_dict = frag.fragment_annotation(output_name_prefix + identifications_file.name,
                                              output_name_prefix + spectrum_file.name,
                                              tolerance,
                                              fragment_types,
                                              charges,
                                              losses,
                                              file_format,
                                              deisotope,
                                              write_file = False)

    # remove written files
    try:
        os.remove(output_name_prefix + spectrum_file.name)
    except Exception as e:
        print("Could not remove file" + output_name_prefix + spectrum_file.name)
    try:
        os.remove(output_name_prefix + identifications_file.name)
    except Exception as e:
        print("Could not remove file" + output_name_prefix + identifications_file.name)

    return fragannot_dict
