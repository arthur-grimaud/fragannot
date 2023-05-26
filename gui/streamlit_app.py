#!/usr/bin/env python3

"""
#####################################################
##                                                 ##
##            -- STREAMLIT MAIN APP --             ##
##                                                 ##
#####################################################
"""

import io
import json
import numpy as np
import pandas as pd

import streamlit as st

from util.converter import JSONConverter
from util.fragannot_call import fragannot_call
from util.redirect import *

@st.cache_data
def dataframe_to_csv_stream(dataframe: pd.DataFrame):
    text = dataframe.to_csv(index = False).encode("utf-8")

    return text

@st.cache_data
def dataframe_to_xlsx_stream(dataframe: pd.DataFrame, sheet_name: str):
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine = "openpyxl") as writer:
        dataframe.to_excel(writer, sheet_name = sheet_name)
        writer.close()

    return buffer

# main page content
def main_page():

    title = st.title("Fragannot - Fragment Ion Annotator")

    general_description = \
    """
    Fragannot is a simple tool implemented in python that annotates fragment ions in MS2 mass spectra. Fragannot returns fragment annotations from a spectra file (in .mgf format) and a corresponding identification file (in .mzid format).
    """
    description = st.markdown(general_description)

    header = st.subheader("Data Import")

    spectrum_file = st.file_uploader("Upload spectrum file:",
                                     type = ["mzml", "mgf"],
                                     help = "Upload a spectrum file to be analyzed in .mzml or .mgf format."
                                    )

    identifications_file = st.file_uploader("Upload identification file:",
                                            type = ["mzid"],
                                            help = "Upload a identification file that contains PSMs of the spectrum file in mzID format."
                                           )

    tolerance = st.number_input("Tolerance in Da:",
                                value = 0.02,
                                help = "Fragment mass tolerance in Dalton.")

    deisotope = st.checkbox("Deisotope spectra",
                            value = True,
                            help = "Deisotope uploaded spectra or not.")

    fions_text = st.markdown("Select which ion types your applied fragmentation method produced:")

    # START Fragment types - List of booleans
    fragannot_ions_col1, fragannot_ions_col2 = st.columns(2)

    with fragannot_ions_col2:

        fions_checkbox_nterm = st.markdown("**N-terminal ions:**")

        fA_ion = st.checkbox("A ions", key = "fA_ion")
        fB_ion = st.checkbox("B ions", value = True, key = "fB_ion")
        fC_ion = st.checkbox("C ions", key = "fC_ion")
        fCdot_ion = st.checkbox("Cdot ions", key = "fCdot_ion")
        fCm1_ion = st.checkbox("C-1 ions", key = "fCm1_ion")
        fCp1_ion = st.checkbox("C+1 ions", key = "fCp1_ion")

    with fragannot_ions_col1:

        fions_checkbox_cterm = st.markdown("**C-terminal ions:**")

        fX_ion = st.checkbox("X ions", key = "fX_ion")
        fY_ion = st.checkbox("Y ions", value = True, key = "fY_ion")
        #fZ_ion = st.checkbox("Z ions", key = "fZ_ion")
        fZdot_ion = st.checkbox("Zdot ions", key = "fZdot_ion")
        fZp1_ion = st.checkbox("Z+1 ions", key = "fZp1_ion")
        fZp2_ion = st.checkbox("Z+2 ions", key = "fZp2_ion")
        fZp3_ion = st.checkbox("Z+3 ions", key = "fZp3_ion")

    fragannot_ion_selection = [fA_ion, fB_ion, fC_ion, fCdot_ion, fCm1_ion, fCp1_ion, fX_ion, fY_ion, fZdot_ion, fZp1_ion, fZp2_ion, fZp3_ion]
    fragannot_ion_names = ["a", "b", "c", "cdot", "c-1", "c+1", "x", "y", "zdot", "z+1", "z+2", "z+3"]
    fragannot_call_ion_selection = []
    for i, sel in enumerate(fragannot_ion_selection):
        if sel:
            fragannot_call_ion_selection.append(fragannot_ion_names[i])

    # END Fragment types

    charges_str = st.text_input("Charges to consider [comma delimited]:",
                                value = "-1, +1",
                                help = "The charges to consider for fragment ions. Multiple entries should be delimited by commas!")
    charges = [charge.strip() for charge in charges_str.split(",")]

    losses_str = st.text_input("Neutral losses to consider [comma delimited]",
                               value = "H2O",
                               help = "Neutral losses to consider for fragment ions. Multiple entries should be delimited by commas!")
    losses = [loss.strip() for loss in losses_str.split(",")]

    l1, l2, center_button, r1, r2 = st.columns(5)

    with center_button:
        run_analysis = st.button("Analyze files!")

    if run_analysis:
        if spectrum_file is not None and identifications_file is not None:
            with st.spinner("Fragannot is running..."):
                with st.expander("Show logging info:"):
                    #with st_stdout("info"):
                    result = fragannot_call(spectrum_file,
                                            identifications_file,
                                            float(tolerance),
                                            fragannot_call_ion_selection,
                                            charges,
                                            losses,
                                            deisotope)
                    converter = JSONConverter()
                    st.session_state["result"] = result
                    st.session_state["dataframes"] = converter.to_dataframes(data = result)
                    status_1 = 0
            if status_1 == 0:
                res_status_1 = st.success("Fragannot finished successfully!")
            else:
                res_status_1 = st.error("Fragannot stopped prematurely! See log for more information!")
        else:
            res_status_1 = st.error("You need to specify a spectrum AND identifications file!")

    if "dataframes" in st.session_state:
        results_header = st.subheader("Download Results")
        csv_1 = st.download_button(label = "Download Fragment-centric CSV!",
                                   data = dataframe_to_csv_stream(st.session_state["dataframes"][0]),
                                   file_name = "fragment_centric.csv",
                                   mime = "text/csv",
                                   help = "Download fragment-centric Fragannot results in CSV format."
                                  )
        csv_2 = st.download_button(label = "Download Spectrum-centric CSV!",
                                   data = dataframe_to_csv_stream(st.session_state["dataframes"][1]),
                                   file_name = "spectrum_centric.csv",
                                   mime = "text/csv",
                                   help = "Download spectrum-centric Fragannot results in CSV format."
                                  )
        #xlsx_1 = st.download_button(label = "Download Fragment-centric XLSX!",
        #                            data = dataframe_to_xlsx_stream(st.session_state["dataframes"][0], "fragment"),
        #                            file_name = "fragment_centric.xlsx",
        #                            mime = "application/vnd.ms-excel",
        #                            help = "Download fragment-centric Fragannot results in XLSX format."
        #                           )
        #xlsx_2 = st.download_button(label = "Download Spectrum-centric XLSX!",
        #                            data = dataframe_to_xlsx_stream(st.session_state["dataframes"][1], "spectrum"),
        #                            file_name = "spectrum_centric.xlsx",
        #                            mime = "application/vnd.ms-excel",
        #                            help = "Download spectrum-centric Fragannot results in XLSX format."
        #                           )
        if "result" in st.session_state:
            json_output = st.download_button(label = "Download JSON!",
                                             data = json.dumps(st.session_state["result"]),
                                             file_name = "result.json",
                                             mime = "text/json",
                                             help = "Download raw Fragannot results in JSON file format."
                                            )

# side bar and main page loader
def main():

    about_str = \
    """
    Fragannot is a simple tool implemented in python that annotates fragment ions in MS2 mass spectra.
    """

    st.set_page_config(page_title = "Fragannot",
                       page_icon = ":test_tube:",
                       layout = "wide",
                       initial_sidebar_state = "expanded",
                       menu_items = {"Get Help": "https://github.com/arthur-grimaud/fragannot/discussions",
                                     "Report a bug": "https://github.com/arthur-grimaud/fragannot/issues",
                                     "About": about_str}
                       )

    title = st.sidebar.title("Fragannot")

    logo = st.sidebar.image("img/logo.gif", caption = "Logo of Fragannot - A golden retriever fetching a fragment ion.")

    doc = st.sidebar.markdown(about_str)

    contact_str = "**Contact:** [Arthur Grimaud](mailto:agrimaud@bmb.sdu.dk), [Veit Schwämmle](veits@bmb.sdu.dk)"
    contact = st.sidebar.markdown(contact_str)

    license_str = "**License:** [MIT License](https://github.com/arthur-grimaud/fragannot/blob/master/LICENSE.md)"
    license = st.sidebar.markdown(license_str)

    project_str = "**Project Page:** [GitHub](https://github.com/arthur-grimaud/fragannot/)"
    project = st.sidebar.markdown(project_str)

    main_page()

if __name__ == "__main__":
    main()
