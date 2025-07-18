import json
import logging
from pathlib import Path

import streamlit as st
from PIL import Image
from streamlit import session_state as state

from background_gene_set import BackgroundGeneSet
from src.enrichment import Enrichment
from src.gene_set import GeneSet
from src.gene_set_library import GeneSetLibrary
from src.ui.helpers import input_example, update_text_widgets
from src.ui.processing import collect_results
from src.ui.rendering import render_results, render_validation
from src.ui.utils import download_link, update_aliases

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)


def main() -> None:
    """
    The main function of the Streamlit application.

    This function sets up the Streamlit app layout, handles user inputs, controls application flow,
    and displays results. It uses various Streamlit components (e.g., st.button, st.dataframe) to
    interact with the user and present information.
    """
    logger.info("Starting the Streamlit app")
    st.sidebar.image(
        Image.open(ROOT / "src" / "static" / "logo.png"),
        caption="Crystal Clear Enrichment Analysis",
    )
    st.sidebar.title("Enrichment analysis")
    st.sidebar.write(
        """This Streamlit app allows users to submit a list of genes and perform enrichment analysis using Gene Ontology pathways. The app displays enriched pathways and GO terms for the submitted genes, along with relevant statistics such as p-values and FDR corrections."""
    )

    st.subheader("Enrichment analysis")

    state.enrich = dict()
    if "results_ready" not in state:
        state.results_ready = False

    state.lib_mapper = update_aliases("libraries")
    state.bg_mapper = update_aliases("backgrounds")
    state.advanced_settings_changed = False
    state.bt_submit_disabled = True

    analysis, advanced_settings = st.tabs(["Analysis", "Advanced settings"])

    with analysis:
        input_gene_set, settings = st.columns([5, 7])
        with input_gene_set:
            st.text_area(
                label="Input a set of genes",
                key="gene_set_input",
                height=438,
                placeholder="Input a gene set",
                label_visibility="collapsed",
            )
            st.text_input(
                label="Gene set name",
                key="gene_set_name",
                placeholder="Input a gene set name",
                label_visibility="collapsed",
            )

            extensions = [".txt"]
            gene_set_files = [
                str(file).replace(f"{ROOT}/data/gene_lists/", "")
                for ext in extensions
                for file in (ROOT / "data" / "gene_lists").rglob(f"*{ext}")
            ]

            capsule_file_selected = st.selectbox(
                "Or select a file from the `data` folder",
                ["Select ..."] + gene_set_files,
                index=0,
                on_change=update_text_widgets,
                key="selected_file",
            )

    with settings:
        state.background_set = st.selectbox(
            "Background gene set", state.bg_mapper.keys()
        )
        st.caption(
            "Specifies the background set of genes. This set validates the input gene set against the chosen organism's genes and serves as a reference for p-value calculations."
        )

        state.libraries = st.multiselect(
            "Select libraries", state.lib_mapper.keys(), default=None
        )

        if ("libraries" in state) and ("lib_mapper" in state):
            state.gene_set_libraries = [
                GeneSetLibrary(
                    str(ROOT / "data" / "libraries" / state.lib_mapper[library]),
                    name=library,
                )
                for library in state.libraries
            ]

        if ("background_set" in state) and ("bg_mapper" in state):
            state.background_gene_set = BackgroundGeneSet(
                str(
                    ROOT
                    / "data"
                    / "backgrounds"
                    / state.bg_mapper[state.background_set]
                )
            )
            if "gene_set_input" in state:
                if state.gene_set_input:
                    state.gene_set = GeneSet(
                        state.gene_set_input.split(),
                        state.background_gene_set.genes,
                        state.gene_set_name,
                    )

    submit, example, placeholder = st.columns([9, 8, 29])
    with submit:
        state.bt_submit_disabled = not all(
            key in state and state[key]
            for key in ["gene_set", "background_gene_set", "gene_set_libraries"]
        )
        bt_submit = st.button("Validate and submit", disabled=state.bt_submit_disabled)

    with example:
        st.button("Input an example", on_click=input_example)

    with advanced_settings:
        n_results = st.slider(
            "Number of results to display", min_value=1, max_value=100, value=10, step=1
        )
        state.p_val_method = st.selectbox(
            "P-value calculation method",
            options=["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
        )
        if state.p_val_method != "Fisher's Exact Test":
            state.advanced_settings_changed = True

        state.bg_custom = st.file_uploader(
            "Upload your background gene set", type=[".txt"]
        )
        if state.bg_custom is not None:
            bg_file = (ROOT / "data" / "backgrounds" / state.bg_custom.name).open("wb")
            bg_file.write(state.bg_custom.getvalue())
            state.advanced_settings_changed = True

        state.libs_custom = st.file_uploader(
            "Upload gene set libraries",
            type=[".gmt"],
            accept_multiple_files=True,
            on_change=update_aliases,
            args=("libraries",),
        )
        for lib_custom in state.libs_custom:
            lib_file = (ROOT / "data" / "libraries" / lib_custom.name).open("wb")
            lib_file.write(lib_custom.getvalue())
            state.advanced_settings_changed = True

        if state.advanced_settings_changed:
            if st.button("Apply settings"):
                logger.info("Running with custom settings")
                st.success("Settings applied")
        else:
            with st.empty():
                st.button("Apply settings", disabled=True)

    if bt_submit:
        logger.info("Validating and submitting genes for enrichment analysis")
        render_validation()
        if state.gene_set_input:
            n_genes = len(state.gene_set_input.split("\n"))
            if (n_genes <= 100) or (n_genes >= 5000):
                if n_genes <= 100:
                    n_warn = "small"
                elif n_genes >= 2000:
                    n_warn = "big"
                s = "s" if str(n_genes)[-1] != 1 else ""
                logger.warning(
                    "You've entered {n_genes} gene{s}, which may be {n_warn} and could affect result accuracy."
                )
                st.warning(
                    f"""You've entered {n_genes} gene{s}, which may be {n_warn} and could affect result accuracy. Consider adjusting p-value or log2 Fold Change.  
    Estimates for the number of DEGs based on comparison type:
    - Similar Conditions (e.g., same cell type, small treatment variations): Dozens to hundreds of DEGs.
    - Moderately Different Conditions (e.g., different cell types, moderate drug treatment): Hundreds to thousands.
    - Highly Different Conditions (e.g., healthy vs. cancerous tissue): Several thousand DEGs."""
                )
            with st.spinner("Calculating enrichment"):
                for gene_set_library in state.gene_set_libraries:
                    logger.info(
                        f"Calculating enrichment results for {gene_set_library.name}"
                    )
                    enrich = Enrichment(
                        state.gene_set,
                        gene_set_library,
                        state.background_gene_set,
                        state.p_val_method,
                    )
                    state.enrich[gene_set_library.name] = enrich
                    with (ROOT / "results" / f"{enrich.name}.json").open(
                        "w"
                    ) as results_snapshot:
                        logger.info(f"Saving {enrich.name}.json")
                        json.dump(enrich.to_snapshot(), results_snapshot)
                logger.info(f"Enrichment results for {gene_set_library.name} are ready")
                state.results_ready = True
        else:
            if not state.gene_set_input:
                logger.error("Please input a newline separated set of genes")
                st.error("Please input a newline separated set of genes")
            if not state.gene_set_libraries:
                logger.error("No libraries were selected for the analysis")
                st.error("No libraries were selected for the analysis")
            if not state.background_gene_set:
                logger.error("No background gene set was selected for the analysis")
                st.error("No background gene set was selected for the analysis")

    if state.results_ready:
        logger.info("Displaying enrichment results")
        st.markdown(
            f'Download all results as {download_link(collect_results(state.enrich), "results", "tsv")}',
            unsafe_allow_html=True,
        )
        for library_name in state.enrich.keys():
            render_results(state.enrich[library_name], library_name, n_results)
        state.results_ready = False

    logger.info("Finishing the Streamlit app")
    return


if __name__ == "__main__":
    main()
