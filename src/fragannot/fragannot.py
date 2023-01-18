# Dependencies
import argparse
from argparse import ArgumentError
from logging import warning
import warnings
import json
import unimod_mapper
import numpy as np
import re
from typing import List  # type hinting
from pyteomics import mass

# Local
from ._version import __version__
from . import constant
from .parser import Parser as Parser

# obj init
um = unimod_mapper.UnimodMapper()

# function called on start
def main(parser=argparse.ArgumentParser()):

    print("-= Starting fragannot =-")

    parser.add_argument("-v", "--version", action="store_true", help="Shows the app version.")
    parser.add_argument(
        "-i",
        "--identification",
        type=str,
        required=True,
        help="Path to spectra identification file in 'mzid' format",
    )
    parser.add_argument(
        "-s", "--spectra", type=str, required=True, help="Path to spectra file in 'mzml' or 'mgf' format"
    )

    parser.add_argument(
        "-t", "--tolerance", type=float, default=0.1, required=False, help="MS2 tolerance in Dalton"
    )

    parser.add_argument(
        "-f",
        "--fragment_types",
        type=str,
        default=["b", "y"],
        nargs="+",
        required=False,
        help="list of fragment ion types to be considered",
    )

    parser.add_argument(
        "-c",
        "--charges",
        type=str,
        default=["+1"],
        nargs="+",
        required=False,
        help="list of fragment charges to be considered",
    )

    parser.add_argument(
        "-l",
        "--losses",
        type=str,
        default=[""],
        nargs="+",
        required=False,
        help="list molecular formula to be considered a possible neutral loss (e.g H20 for water loss)",
    )

    args = parser.parse_args()

    if args.version:
        return __version__
    else:

        print_parameters(vars(args))

        fragment_annotation(
            ident_file=args.identification,
            spectra_file=args.spectra,
            tolerance=args.tolerance,
            fragment_types=args.fragment_types,
            charges=args.charges,
            losses=args.losses,
        )


def print_parameters(args):

    for arg, val in args.items():
        print(f"{arg} : {val} \n")


def fragment_annotation(ident_file, spectra_file, tolerance, fragment_types, charges, losses):

    P = Parser()

    psms = P.read(spectra_file, ident_file)
    i = 0

    psms_json = []

    for psm in psms:

        print(i)
        theoretical_fragment_code = compute_theoretical_fragments2(
            sequence_length=len(psm.peptidoform.sequence),
            fragment_types=fragment_types,
            charges=charges,
            neutral_losses=losses,
        )

        theoretical_fragment_dict = {
            f: theoretical_mass_to_charge(f, psm.peptidoform) for f in theoretical_fragment_code
        }

        annotation_mz, annotation_code = match_fragments(
            psm.spectrum["mz"], theoretical_fragment_dict, tolerance=tolerance
        )

        psm.spectrum["intensity"] = psm.spectrum["intensity"].tolist()
        psm.spectrum["mz"] = psm.spectrum["mz"].tolist()
        psm.spectrum["theoretical_mz"] = annotation_mz
        psm.spectrum["theoretical_code"] = annotation_code

        # add to json

        psms_json.append(
            {
                "sequence": psm.peptidoform.sequence,
                "proforma": psm.peptidoform.proforma,
                "annotation": psm.spectrum,
            }
        )

    with open(P.output_fname, "w", encoding="utf8") as f:
        json.dump(psms_json, f)


# Function


def compute_theoretical_fragments2(
    sequence_length: int,
    fragment_types: List[str],
    charges: List[int] = [1],
    neutral_losses: List[str] = [],
    internal: bool = True,
) -> List[str]:

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
        internal = [n_term_ion + ":" + c_term_ion for n_term_ion in n_term_ions for c_term_ion in c_term_ions]
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
            internal_frag + nl for internal_frag in internal_frags_with_charges for nl in neutral_losses_str
        ]

    return n_term_frags_with_nl + c_term_frags_with_nl + internal_frags_with_nl


def theoretical_mass_to_charge(fragment_code, peptidoform):

    start, end, ion_cap_start, ion_cap_end, charge, formula = parse_fragment_code(fragment_code)

    # peptide and modification mass
    sequence = []
    mods = []
    for aa, mod in peptidoform.parsed_sequence[start - 1 : end]:
        sequence.append(aa)
        if not mod is None:
            mods.extend([m.mass for m in mod])

    # mass AA sequence
    P = mass.fast_mass(sequence="".join(sequence))
    # mass modifications
    M = sum(mods)
    # mass start ion cap
    SI = constant.ion_cap_delta_mass[ion_cap_start]
    # mass end ion cap
    EI = constant.ion_cap_delta_mass[ion_cap_end]
    # hydrogen mass
    H = mass.calculate_mass("H", absolute=True)
    # loss mass
    L = mass.calculate_mass(formula, absolute=True)

    # Calculate fragment mass
    fragment_mass = (P + M + SI + EI + (H * charge) - L) / np.abs(charge)

    return fragment_mass


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


def match_fragments(fragment_mz, frag_theo, tolerance):

    fragment_theoretical_mz = []
    fragment_theoretical_code = []

    for mz in fragment_mz:
        closest_frag = min(frag_theo.items(), key=lambda x: abs(x[1] - mz))
        if abs(mz - closest_frag[1]) <= tolerance:
            fragment_theoretical_mz.append(closest_frag[1])
            fragment_theoretical_code.append(closest_frag[0])
        else:
            fragment_theoretical_mz.append(None)
            fragment_theoretical_code.append(None)

    return (fragment_theoretical_mz, fragment_theoretical_code)
