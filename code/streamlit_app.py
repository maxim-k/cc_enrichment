import base64
import json
import logging
from math import log10
from typing import Dict
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st
from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from PIL import Image
from streamlit import session_state as state

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)


def update_aliases(directory: str, alias_file: str = "alias.json") -> Dict[str, str]:
    """
    Update the aliases file with the current set of files in the specified directory.

    This function scans a specified directory for files, keeping an 'alias.json' file that contains
    a mapping of simplified names to actual file names. This allows for easier reference to the files
    in that directory. If new files are added to the directory, they are added to the 'alias.json'. If files
    are missing, they are removed from 'alias.json'.

    :param directory: The directory to scan for files.
    :param alias_file: The name of the alias file, defaults to 'alias.json'
    :return: A dictionary containing the aliases, with keys being the simplified names and values being the actual file names.
    """
    logger.info(f"Updating aliases for directory: {directory}")
    aliases_path = ROOT / "data" / directory / alias_file

    if Path(aliases_path).is_file():
        try:
            with open(aliases_path, "r") as file:
                alias = json.load(file)
        except (FileNotFoundError, json.JSONDecodeError):
            logger.warning(f"Failed to load aliases from {aliases_path}")
            st.warning(f"Failed to load aliases from {aliases_path}")

    files = [
        f
        for f in (ROOT / "data" / directory).iterdir()
        if f.is_file()
    ]

    # Remove 'alias.json' from the list of files
    if Path(aliases_path) in files:
        files.remove(Path(aliases_path))

    # Add a record if a file is not in aliases
    for file in files:
        if file.name not in alias.values():
            alias[file.stem] = file.name

    # Delete a record from aliases if there's no corresponding file
    aliases_keys_to_delete = [key for key in alias if alias[key] not in [file.name for file in files]]

    for key in aliases_keys_to_delete:
        del alias[key]

    with open(aliases_path, "w") as file:
        json.dump(alias, file, indent=4)

    return alias


def render_table(result: pd.DataFrame) -> None:
    """
    Render a styled DataFrame within the Streamlit app.

    This function takes a DataFrame containing results data and applies custom formatting before
    displaying it within the Streamlit app using the `st.dataframe` function. Custom formatting includes
    specific number formatting for 'p-value' and 'fdr' columns.

    :param result: The DataFrame containing results data to display.
    """

    logger.info("Rendering DataFrame in Streamlit app.")
    def custom_format(n):
        if n > 0.001:
            return f"{n:.3f}"
        else:
            return f"{n:.3e}"

    st.dataframe(
        result.style.format({"p-value": custom_format, "fdr": custom_format}),
        use_container_width=True,
        column_config={
            "rank": None,
            "term": "Term",
            "overlap_size": "Overlap size",
            "overlap": "Overlap",
            "description": None,
            "p-value": "P-value",
            "fdr": "FDR",
        },
    )


def render_barchart(result: pd.DataFrame) -> None:
    """
    Render a bar chart visualization of the results data within the Streamlit app.

    This function takes a DataFrame containing results data, specifically terms and p-values, and creates
    a bar chart using Plotly Express, which is then displayed in the Streamlit app using the `st.plotly_chart` function.

    :param result: The DataFrame containing results data to visualize.
    """
    logger.info("Rendering bar chart in Streamlit app.")
    bar = result[["term", "p-value"]]
    bar.loc[:, "p-value"] = bar.loc[:, "p-value"].apply(lambda x: -1 * log10(x))
    bar = bar.sort_values(by=["p-value"])
    bar.columns = ["term", "-log10(p-value)"]
    fig = px.bar(
        bar,
        x="-log10(p-value)",
        y="term",
        orientation="h",
        labels={"-log10(p-value)": "−log₁₀(p‐value)", "term": "Term"},
    )
    st.plotly_chart(fig, use_container_width=True)


def download_link(val: str, filename: str, extension: str) -> str:
    """
    Create a download link for a file with the given content, filename, and extension.

    This function generates a download link that, when clicked, will download the file with the given
    content, filename, and extension. The content is encoded in base64 and the link is created using an
    HTML 'a' tag with a 'download' attribute.

    :param val: The content of the file to be downloaded.
    :param filename: The name of the file, without the extension.
    :param extension: The file extension (e.g., 'tsv', 'json').
    :return: An HTML string containing the download link.
    """
    logger.info(f"Creating download link for file: {filename}.{extension}")
    b64 = base64.b64encode(val.encode("utf-8"))
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{filename}.{extension}">{extension}</a>'


def collect_results(results: Dict) -> str:
    """
    Concatenate all enrichment results. For each result adds a column with a library name.

    :param results: The dictionary containing enrichment results.
    """
    logger.info("Concatenating all enrichment results.")
    results_concat = []
    for library_name in results.keys():
        result = results[library_name].to_dataframe()
        result.insert(0, 'Library', library_name)
        results_concat.append(result)

    return pd.concat(results_concat, ignore_index=True).to_csv(sep="\t", index=False)


def render_results(result: Enrichment, file_name: str, n_results: int = 10) -> None:
    """
    Render a results section within the Streamlit app.

    This function processes and visualizes the result data within the Streamlit app. It provides
    a table visualization and a bar chart visualization of the top results. Additionally, it
    allows the user to download the results in various formats.

    :param result: The DataFrame containing results data to display.
    :param file_name: The name of the file to be used for downloading results.
    :param n_results: Numbers of results to display
    """
    logger.info(f"Rendering results for file: {file_name}")
    result_df = result.to_dataframe().head(n_results)
    result_df = result_df.set_index("rank")
    table, bar = st.tabs(["Results", "Bar chart"])
    with table:
        render_table(result_df)
    with bar:
        render_barchart(result_df)

    st.markdown(
        f'Download results as {download_link(result.to_tsv(), file_name, "tsv")}, {download_link(result.to_json(), file_name, "json")}',
        unsafe_allow_html=True,
    )


def render_validation() -> None:
    """
    Validate and render the gene set information.

    This function checks the `gene_set` in the session state for duplicates and invalid entries.
    It then provides a feedback to the user in the Streamlit app on the validation results.
    """
    logger.info("Validating and rendering the gene set information.")
    if "gene_set" in state:
        total = state.gene_set.size
        dups = len(state.gene_set.validation["duplicates"])
        non_gene = len(state.gene_set.validation["non_hgnc"])
        if dups:
            dups_st = f", ⚠️ {dups} duplicates"
        else:
            dups_st = ""
        if non_gene:
            non_gene_st = f", ⛔ {non_gene} non valid"
        else:
            non_gene_st = ""
        caption = f"{total} genes{dups_st}{non_gene_st}"
        with st.expander(caption):
            if dups:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)["duplicates"],
                    column_config={
                        "duplicates": st.column_config.ListColumn(
                            dups_st[2:], width="large"
                        ),
                    },
                    hide_index=True,
                )
            if non_gene:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)["non_hgnc"],
                    column_config={
                        "non_hgnc": st.column_config.ListColumn(
                            non_gene_st[2:], width="large"
                        ),
                    },
                    hide_index=True,
                )
        state.bt_submit_disabled = False


def input_example() -> None:
    """
    Set the example input for the Streamlit app.

    This function loads the example gene set and populates the `gene_set_input` and
    `gene_set_name` in the session state with the content and name of the example respectively.
    """
    logger.info("Setting the example input for the Streamlit app.")
    # Callback because that's the only way it works
    state.gene_set_input = (ROOT / "code" / "static" / "example_gene_list.txt").read_text()
    state.gene_set_name = "Example gene set"


def main() -> None:
    """
    The main function of the Streamlit application.

    This function sets up the Streamlit app layout, handles user inputs, controls application flow,
    and displays results. It uses various Streamlit components (e.g., st.button, st.dataframe) to
    interact with the user and present information.
    """
    logger.info("Starting the Streamlit app")
    st.sidebar.image(
        Image.open(f"{ROOT}/code/static/CO_logo_135x72.png"), caption="Code Ocean"
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

        with settings:
            state.background_set = st.selectbox("Background gene set", state.bg_mapper.keys())
            st.caption(
                "Specifies the background set of genes. This set validates the input gene set against the chosen organism's genes and serves as a reference for p-value calculations."
            )

            state.libraries = st.multiselect(
                "Select libraries",
                state.lib_mapper.keys(),
                default=None,
            )

            if ("libraries" in state) and ("lib_mapper" in state):
                state.gene_set_libraries = [
                    GeneSetLibrary(
                        str(ROOT / "data" / "libraries" / state.lib_mapper[library]), name=library
                    )
                    for library in state.libraries
                ]

            if ("background_set" in state) and ("bg_mapper" in state):
                state.background_gene_set = BackgroundGeneSet(
                    str(ROOT / "data" / "backgrounds" / state.bg_mapper[state.background_set])
                )
                if "gene_set_input" in state:
                    state.bt_submit_disabled = False
                    if state.gene_set_input:
                        state.gene_set = GeneSet(
                            state.gene_set_input.split(), state.background_gene_set.genes, state.gene_set_name
                        )

        submit, example, placeholder = st.columns([8, 8, 30])
        with submit:
            bt_submit = st.button(
                "Validate and submit", disabled=state.bt_submit_disabled
            )

        with example:
            st.button("Input an example", on_click=input_example)
    with advanced_settings:
        n_results = st.slider(
            "Number of results to display", min_value=1, max_value=100, value=10, step=1
        )
        p_val_method = st.selectbox(
            "P-value calculation method",
            options=["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
        )
        bg_custom = st.file_uploader("Upload your background gene set", type=[".txt"])
        if bg_custom is not None:
            bg_file = (ROOT / "data" / "backgrounds" / bg_custom.name).open("wb")
            bg_file.write(bg_custom.getvalue())
            state.advanced_settings_changed = True

        libs_custom = st.file_uploader(
            "Upload gene set libraries", type=[".gmt"], accept_multiple_files=True, on_change=update_aliases, args=("libraries", )
        )
        for lib_custom in libs_custom:
            lib_file = (ROOT / "data" / "libraries" / lib_custom.name).open("wb")
            lib_file.write(lib_custom.getvalue())
            state.advanced_settings_changed = True

        if state.advanced_settings_changed:
            if st.button("Apply settings"):
                st.success("Settings applied")
        else:
            with st.empty():
                st.button("Apply settings", disabled=True)

    if bt_submit:
        logger.info("Validating and submitting genes for enrichment analysis")
        render_validation()
        if state.gene_set_input:
            n_genes = len(state.gene_set_input.split("\n"))
            if (n_genes <= 10) or (n_genes >= 5000):
                if n_genes <= 100:
                    n_warn = "small"
                elif n_genes >= 2000:
                    n_warn = "big"
                s = "s" if str(n_genes)[-1] != 1 else ""
                st.warning(
                    f"""You've entered {n_genes} gene{s}, which may be {n_warn} and could affect result accuracy. Consider adjusting p-value or log2 Fold Change.  
Estimates for the number of DEGs based on comparison type:
- Similar Conditions (e.g., same cell type, small treatment variations): Dozens to hundreds of DEGs.
- Moderately Different Conditions (e.g., different cell types, moderate drug treatment): Hundreds to thousands.
- Highly Different Conditions (e.g., healthy vs. cancerous tissue): Several thousand DEGs."""
                )
            with st.spinner("Calculating enrichment"):
                logger.info("Calculating enrichment for the submitted genes")
                for gene_set_library in state.gene_set_libraries:
                    enrich = Enrichment(
                        state.gene_set,
                        gene_set_library,
                        state.background_gene_set,
                        p_val_method,
                    )
                    state.enrich[gene_set_library.name] = enrich
                    with (ROOT / "results" / f"{enrich.name}.json").open("w") as results_snapshot:
                        json.dump(enrich.to_snapshot(), results_snapshot)
                state.results_ready = True
        else:
            if not state.gene_set_input:
                st.error("Please input a newline separated set of genes")
            if not state.gene_set_libraries:
                st.error("No libraries were selected for the analysis")

    if state.results_ready:
        logger.info("Displaying enrichment results")
        st.divider()
        st.markdown(
            f'Download all results as {download_link(collect_results(state.enrich), "results", "tsv")}',
            unsafe_allow_html=True,
        )
        for library_name in state.enrich.keys():
            st.subheader(library_name)
            render_results(state.enrich[library_name], library_name, n_results)
    logger.info("Ending the Streamlit app")
    return


if __name__ == "__main__":
    main()
