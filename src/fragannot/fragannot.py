import argparse
from ._version import __version__

# Generic
from argparse import ArgumentError
from logging import warning
import warnings
import json

# Development
from pprintpp import pprint

# MS file parsers
from pyteomics import mzid
from pyteomics import mgf
from pyteomics import mzml
from pyteomics import mass

#
import unimod_mapper
import re
import numpy as np

#
from . import constant
from .parser import Parser as Parser
import spectrum_utils.spectrum as sus

# type hinting
from typing import List

#
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

    args = parser.parse_args()

    if args.version:
        return __version__
    else:

        fragment_annotation(ident_file=args.identification, spectra_file=args.spectra)
        return "end return"


def fragment_annotation(ident_file, spectra_file):

    P = Parser()

    psms = P.read(spectra_file, ident_file)
    i = 0

    psms_json = []

    for psm in psms:

        print(i)
        theoretical_fragment_code = compute_theoretical_fragments2(
            sequence_length=len(psm.peptidoform.sequence), fragment_types=["b", "y"]
        )

        theoretical_fragment_dict = {
            f: theoretical_mass_to_charge(f, psm.peptidoform) for f in theoretical_fragment_code
        }

        annotation_mz, annotation_code = match_fragments(
            psm.spectrum["mz"], theoretical_fragment_dict, tolerance=0.1
        )

        psm.spectrum["intensity"] = psm.spectrum["intensity"].tolist()
        psm.spectrum["mz"] = psm.spectrum["mz"].tolist()
        psm.spectrum["theoretical_mz"] = annotation_mz
        psm.spectrum["theoretical_code"] = annotation_code

        # add to json
        print(psm.spectrum)

        psms_json.append(
            {
                "sequence": psm.peptidoform.sequence,
                "proforma": psm.peptidoform.proforma,
                "annotation": psm.spectrum,
            }
        )

        if i == 1000:
            break
        i += 1

    with open("data.json", "w", encoding="utf8") as f:
        json.dump(psms_json, f)


# READERS/PARSERS #


def get_modification_proforma(sequence, modification, mod_mass_to_name=None):
    """Returns the proteoform modification in the proforma format"""
    with warnings.catch_warnings():  # catch warnings from unimod
        if modification == []:
            return sequence
        else:
            sequence_list = list(sequence.strip(" "))

            for mod in reversed(modification):

                modMass = mod["monoisotopicMassDelta"]

                if mod_mass_to_name is None:  # if called without mod to mass dict
                    modName = um.id_to_name(um.mass_to_ids(modMass, decimals=2)[0])[0]
                elif modMass not in mod_mass_to_name:  # if called with mod to mass dict but mass not added
                    modName = um.id_to_name(um.mass_to_ids(modMass, decimals=2)[0])[0]
                    mod_mass_to_name[modMass] = modName  # mass found in mass
                else:
                    modName = mod_mass_to_name[modMass]

                sequence_list.insert(int(mod["location"]), "[{0}]".format(modName))

            return "".join(sequence_list), mod_mass_to_name


def read_mzml(self, spectra_file):
    """
    Add informations from an mzml file to the spectrum objects.

    This method reads the spectra in an mzml file and adds that information to the Spetrum object in self.spectra

    Parameters
    ----------
    spectra_file : str
        The path to the .mgf file to be read.

    Returns
    -------
    """
    self.spectra_file = spectra_file  # Store File name that has been read
    mzml_obj = mzml.read(spectra_file)

    print("---Loading spectrum data from: {0}---".format(self.spectra_file))
    with alive_bar(0) as bar:
        sp_to_del = []
        for spec_id in self.spectra:  # Take into account how DBSEs store spectra ids
            print(spec_id)

            if self.dbse == "mascot":
                index = str(int(spec_id.split("=")[1]) + 1)
            if self.dbse == "comet":
                index = spec_id.split("=")[1]

            try:
                spec_mzml = mzml_obj.get_by_id(index, id_key="index")  # need to be splited
            except (KeyError):
                print("Key error spectrum not found")

            if spec_mzml["ms level"] != 1:
                self.spectra[spec_id].set_spec_data_mzml(spec_mzml)
            else:
                sp_to_del.append(spec_id)
            bar()

            # delete ms1
        for spec_id in sp_to_del:
            del self.spectra[spec_id]

    pass


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

    charges_str = ["(" + str(charge) + ")" if charge < 0 else "(+" + str(charge) + ")" for charge in charges]
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
    #
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


def compute_theoretical_fragments(sequence, modifications, ionTypes):
    """Returns and set a list of m/z of fragment ions  and informations on the type/position of each fragments for a given peptidoform/proteoform"""

    ion_formulas = constant.ion_formulas

    frag_masses = {}

    # fragments masses:
    for ion_type in ionTypes:
        frag_masses_iontype = {}

        if "-" in ion_type:  # Internal Fragment
            # sum of all modification masses present in the internal fragment
            sum_mods = lambda modifications, i, j: sum(
                [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"] <= j]
            )  # sum mods delta for internal fragm ions
            # get all sub string of the peptide sequence
            sub_sequences = [
                (
                    sequence[i - 1 : j],
                    i,
                    j,
                    ion_type,
                    [mod["location"] for mod in modifications if i <= mod["location"] <= j],
                )
                for i in range(1, len(sequence))
                for j in range(i + 1, len(sequence))
            ]
            # compute internal frag masses
            frag_masses_iontype.update(
                {
                    ",".join(str(s) for s in seq[1:4]): round(
                        mass.fast_mass(
                            sequence=seq[0],
                            ion_type=ion_type,
                            ion_comp=ion_formulas[ion_type],
                        )
                        + sum_mods(modifications, seq[1], seq[2]),
                        4,
                    )
                    for seq in sub_sequences
                }
            )

        else:  # Terminal Fragment
            if any(i in ion_type for i in ["a", "b", "c"]):  # Nterm
                sum_mods = lambda modifications, i, j: sum(
                    [mod["monoisotopicMassDelta"] for mod in modifications if mod["location"] <= j]
                )
                sub_sequences = [
                    (
                        sequence[:j],
                        1,
                        j,
                        ion_type,
                        [mod["location"] for mod in modifications if mod["location"] <= j],
                    )
                    for j in range(2, len(sequence))
                ]
                frag_masses_iontype.update(
                    {
                        ",".join(str(s) for s in seq[1:4]): round(
                            mass.fast_mass(
                                sequence=seq[0],
                                ion_type=ion_type,
                                ion_comp=ion_formulas[ion_type],
                            )
                            + sum_mods(modifications, seq[1], seq[2]),
                            4,
                        )
                        for seq in sub_sequences
                    }
                )

            else:  # Cterm
                sum_mods = lambda modifications, i, j: sum(
                    [mod["monoisotopicMassDelta"] for mod in modifications if i <= mod["location"]]
                )
                sub_sequences = [
                    (
                        sequence[i - 1 :],
                        i,
                        len(sequence),
                        ion_type,
                        [mod["location"] for mod in modifications if i <= mod["location"]],
                    )
                    for i in range(1, len(sequence) + 1)
                ]
                frag_masses_iontype.update(
                    {
                        ",".join(str(s) for s in seq[1:4]): round(
                            mass.fast_mass(
                                sequence=seq[0],
                                ion_type=ion_type,
                                ion_comp=ion_formulas[ion_type],
                            )
                            + sum_mods(modifications, seq[1], seq[2]),
                            4,
                        )
                        for seq in sub_sequences
                    }
                )

        frag_masses.update(frag_masses_iontype)

    return frag_masses


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
