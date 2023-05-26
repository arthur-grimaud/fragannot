# Dependencies
import json
import numpy as np
import re
from pyteomics import mass
from itertools import tee
from random import uniform
from pyteomics import parser
import os
import ms_deisotope

# type hinting
from typing import List
from typing import Dict
from typing import Any

from fragannot import constant
from fragannot.parser import Parser as Parser

import codon

class FragannotCodon:
    def __init__(self):
        pass

    def fragment_annotation(
        self,
        ident_file: str,
        spectra_file: str,
        tolerance: float,
        fragment_types: List[str],
        charges: List[str],
        losses: List[str],
        file_format: str,
        deisotope: bool,
        write_file: bool = True) -> List[Dict[str, Any]]:

        return fragment_annotation(ident_file, spectra_file, tolerance,
                                   fragment_types, charges, losses, file_format,
                                   deisotope, write_file)

def fragment_annotation(
    ident_file: str,
    spectra_file: str,
    tolerance: float,
    fragment_types: List[str],
    charges: List[str],
    losses: List[str],
    file_format: str,
    deisotope: bool,
    write_file: bool = True) -> List[Dict[str, Any]]:
    """
    Annotate theoretical and observed fragment ions in a spectra file.

    Parameters:
    ----------
    ident_file : str
        Filename of an identification file
    spectra_file : str
        Filename of a spectra file
    tolerance : float
        Tolerance value in ppm for fragment matching
    fragment_types : list
        List of fragment types (fragment type must be defined in constant.py ion_cap_formula)
    charges : list
        List of charges (e.g ["+1", "-2"])
    losses : list
        List of neutral losses molecular formula (e.g ["H2O"])
    file_format : str
        String indicating the file format of the input files

    Returns:
    -------
    None
    """

    P = Parser()

    psms = P.read(spectra_file, ident_file, file_format = file_format)
    i = 0

    psms_json = []

    for psm in psms:

        if (i + 1) % 100 == 0:
            print(f"{i + 1} spectra annotated")

        if charges == "auto":  # if charges to consider not specified: use precursor charge as max charge
            charges_used = range(1, abs(psm.get_precursor_charge()), 1)
        else:
            charges_used = charges

        theoretical_fragment_code = compute_theoretical_fragments(
            sequence_length = len(psm.peptidoform.sequence),
            fragment_types = fragment_types,
            charges = charges_used,
            neutral_losses = losses,
        )

        theoretical_fragment_dict = {
            f: theoretical_mass_to_charge(f, psm.peptidoform) for f in theoretical_fragment_code
        }

        if deisotope:  # deisotoping TODO check desotoping method to optimize
            mzs, intensities = deisotope_peak_list(
                psm.spectrum["mz"].tolist(), psm.spectrum["intensity"].tolist()
            )
        else:
            mzs = psm.spectrum["mz"]
            intensities = psm.spectrum["intensity"].tolist()

        annotation_mz, annotation_code, annotation_count = match_fragments(
            mzs, theoretical_fragment_dict, tolerance = tolerance
        )

        psm.spectrum["intensity"] = intensities
        psm.spectrum["mz"] = mzs
        psm.spectrum["theoretical_mz"] = annotation_mz
        psm.spectrum["theoretical_code"] = annotation_code
        psm.spectrum["matches_count"] = annotation_count

        psms_json.append(
            {
                "sequence": psm.peptidoform.sequence,
                "proforma": psm.peptidoform.proforma,
                "annotation": psm.spectrum,
                "spectrum_id": psm.spectrum_id,
                "identification_score": psm.score,
                "rank": psm.rank,
                # "precursor_charge": int(psm.get_precursor_charge()),
                "precursor_intensity": 666,
            }
        )
        i += 1

    if write_file:
        with open(P.output_fname, "w", encoding="utf8") as f:
            json.dump(psms_json, f)

    return psms_json

def deisotope_peak_list(mzs: List[float], intensities: List[float]) -> List[List[float]]:
    peaks = ms_deisotope.deconvolution.utils.prepare_peaklist(zip(mzs, intensities))
    deconvoluted_peaks, targeted = ms_deisotope.deconvolute_peaks(
        peaks, averagine = ms_deisotope.peptide, scorer = ms_deisotope.MSDeconVFitter(10.0), verbose = True
    )
    mzs = [p.mz for p in deconvoluted_peaks.peaks]
    intensities = [p.intensity for p in deconvoluted_peaks.peaks]

    return mzs, intensities

@codon.jit(pyvars=["constant"], debug=True)
def compute_theoretical_fragments(
    sequence_length: int,
    fragment_types: List[str],
    charges: List[int] = [-1],
    neutral_losses: List[str] = [],
    internal: bool = True) -> List[str]:

    ion_directions = constant.ion_direction

    n_term_ions = [ion_type for ion_type in fragment_types if ion_directions[ion_type] == "n-term"]
    c_term_ions = [ion_type for ion_type in fragment_types if ion_directions[ion_type] == "c-term"]

    n_term = ["t:" + ion_type for ion_type in n_term_ions]
    c_term = [ion_type + ":t" for ion_type in c_term_ions]

    # terminal fragments
    n_term_frags = [
        n_term_frag + "@1:" + str(i + 1) for n_term_frag in n_term for i in range(sequence_length - 1)
    ]
    c_term_frags = [
        c_term_frag + "@" + str(i) + ":" + str(sequence_length)
        for c_term_frag in c_term
        for i in range(2, sequence_length + 1)
    ]

    charges_str = [
        "(" + str(int(charge)) + ")" if int(charge) < 0 else "(+" + str(int(charge)) + ")"
        for charge in charges
    ]

    n_term_frags_with_charges = [
        n_term_frag + charge for n_term_frag in n_term_frags for charge in charges_str
    ]

    c_term_frags_with_charges = [
        c_term_frag + charge for c_term_frag in c_term_frags for charge in charges_str
    ]

    neutral_losses_str = ["[" + nl + "]" for nl in neutral_losses]
    neutral_losses_str.append("")
    n_term_frags_with_nl = [
        n_term_frag + nl for n_term_frag in n_term_frags_with_charges for nl in neutral_losses_str
    ]
    c_term_frags_with_nl = [
        c_term_frag + nl for c_term_frag in c_term_frags_with_charges for nl in neutral_losses_str
    ]

    internal_frags_with_nl = []

    if internal:
        # internal fragments
        internal = [
            n_term_ion + ":" + c_term_ion for n_term_ion in n_term_ions for c_term_ion in c_term_ions
        ]
        internal_pos = [
            str(i) + ":" + str(j)
            for i in range(2, sequence_length)
            for j in range(2, sequence_length)
            if i <= j
        ]
        internal_frags = [
            internal_ions + "@" + internal_positions
            for internal_ions in internal
            for internal_positions in internal_pos
        ]

        internal_frags_with_charges = [
            internal_frag + charge for internal_frag in internal_frags for charge in charges_str
        ]

        internal_frags_with_nl = [
            internal_frag + nl
            for internal_frag in internal_frags_with_charges
            for nl in neutral_losses_str
        ]

    return n_term_frags_with_nl + c_term_frags_with_nl + internal_frags_with_nl

@codon.jit(pyvars=["parse_fragment_code", "parser", "mass", "constant"])
def theoretical_mass_to_charge(fragment_code: str, peptidoform) -> float:

    start, end, ion_cap_start, ion_cap_end, charge, formula = parse_fragment_code(fragment_code)

    # peptide and modification mass
    sequence = []
    mods = []
    for aa, mod in peptidoform.parsed_sequence[start - 1 : end]:
        sequence.append(aa)
        if not mod is None:
            mods.extend([m.mass for m in mod])

    # mass AA sequence
    ps = parser.parse("".join(sequence), show_unmodified_termini=True)
    P = mass.calculate_mass(parsed_sequence=ps)
    # mass modifications
    M = sum(mods)
    # mass start ion cap
    SI = constant.ion_cap_delta_mass[ion_cap_start]
    # mass end ion cap
    EI = constant.ion_cap_delta_mass[ion_cap_end]
    # hydrogen mass
    H = 1.00784
    # loss mass
    L = mass.calculate_mass(formula, absolute=True)

    # Calculate fragment mass
    fragment_mass = (P + M + SI + EI + (H * charge) - L) / abs(charge)

    return fragment_mass

@codon.jit(pyvars=["re"])
def parse_fragment_code(fragment_code: str):

    # test if fragment code format is valid*
    fragment_code_pattern = re.compile(".+(:).+(@)[0-9]+(:)[0-9]+(\()(\+|\-)[0-9](\))(\[(.*?)\])?")
    if bool(fragment_code_pattern.match(fragment_code)) == False:
        raise RuntimeError("Incorrect fragment code format: {0}".format(fragment_code))

    ## Parse fragment code

    start, end = [
        int(i) for i in re.search("(?<=\@)(.*?)(?=\()", fragment_code).group(1).split(":")
    ]  # Get start and end amino acid indexes
    ion_cap_start, ion_cap_end = [
        str(i) for i in re.search("^(.*?)(?=\@)", fragment_code).group(1).split(":")
    ]  # Get start and end ion caps name
    charge = int(re.search("(?<=\()(.*?)(?=\))", fragment_code).group(1))  # get charge state
    formula = re.search("(?<=\[)(.*?)(?=\])", fragment_code)
    if formula == None:
        formula = ""
    else:
        formula = str(re.search("(?<=\[)(.*?)(?=\])", fragment_code).group(1))

    return start, end, ion_cap_start, ion_cap_end, charge, formula

@codon.jit(pyvars=["re", "tee"])
def match_fragments(exp_mz, theo_frag, tolerance: float):

    theo_frag = [[k, v] for k, v in sorted(theo_frag.items(), key = lambda item: item[1])]

    re_term = re.compile(r"^t:|:t")

    iter_2, last_match = tee(iter(theo_frag))

    d = {}

    fragment_theoretical_code = []
    fragment_theoretical_mz = []
    fragment_theoretical_nmatch = []

    for i in exp_mz:
        d.setdefault(i, [])
        found = False
        while True:
            j = next(iter_2, (None, None))

            # print(j)
            if j[1] is None:
                break
            if abs(i - j[1]) <= tolerance:
                k = [j[0], j[1], abs(i - j[1])]

                d[i].append(k)
                if not found:
                    iter_2, last_match = tee(iter_2)
                    found = True
            else:
                if found:
                    break

        fragment_theoretical_nmatch.append(len(d[i]))
        if len(d[i]) > 0:

            closest = None
            for frag in d[i]:

                if re_term.search(frag[0]):  # Prioritize annotation of terminal ions
                    closest = frag
                    break

            if closest is None:
                closest = min(
                    d[i], key = lambda t: t[2]
                )  # add the only the annotation with the lowest mass error

            fragment_theoretical_code.append(closest[0])
            fragment_theoretical_mz.append(closest[1])
        else:
            fragment_theoretical_code.append(None)
            fragment_theoretical_mz.append(None)

        iter_2, last_match = tee(last_match)

    return (fragment_theoretical_mz, fragment_theoretical_code, fragment_theoretical_nmatch)
