from math import log10
import json
import os
import base64

import plotly.express as px
import streamlit as st
from PIL import Image
import pandas as pd

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from parallel_enrichment import Enrichment
from streamlit import session_state as state

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)

ROOT = "/Users/codeocean/PycharmProjects/co_enrichment/"

UPLOAD_FOLDER = f"{ROOT}results/upload"


def update_aliases(directory, alias_file="alias.json"):
    aliases_path = os.path.join(f"{ROOT}data/{directory}", alias_file)
    aliases = {}

    # Check if aliases.json exists; if so, load its content.
    if os.path.isfile(aliases_path):
        with open(aliases_path, 'r') as file:
            aliases = json.load(file)

    # List all files in the directory.
    files = [f for f in os.listdir(f"{ROOT}data/{directory}") if
             os.path.isfile(os.path.join(f"{ROOT}data/{directory}", f))]

    # Remove 'alias.json' from the list of files.
    if alias_file in files:
        files.remove(alias_file)

    # Check for files that are not in the aliases dictionary and add them.
    for file in files:
        if file not in aliases.values():
            # Extract filename without extension.
            file_name_without_ext = os.path.splitext(file)[0]
            aliases[file_name_without_ext] = file

    # Save the updated aliases back to aliases.json.
    with open(aliases_path, 'w') as file:
        json.dump(aliases, file, indent=4)

    return aliases


lib_mapper = update_aliases("libraries")
bg_mapper = update_aliases("backgrounds")


def render_table(result):
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


def render_barchart(result):
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


def download_link(val, filename, extension):
    b64 = base64.b64encode(val.encode('utf-8'))
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{filename}.{extension}">{extension}</a>'


def render_results(result):
    result_df = result.to_dataframe().head(10)
    result_df = result_df.set_index("rank")
    table, bar = st.tabs(["Results", "Bar chart"])
    with table:
        render_table(result_df)
    with bar:
        render_barchart(result_df)

    st.markdown(
        f'Download results as {download_link(result.to_tsv(), "result", "tsv")}, {download_link(result.to_json(), "result", "json")}',
        unsafe_allow_html=True)


def render_validation():
    if "gene_set" in state:
        total = state.gene_set.size
        dups = len(state.gene_set.validation['duplicates'])
        non_gene = len(state.gene_set.validation['non_hgnc'])
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
                    pd.json_normalize(state.gene_set.validation)['duplicates'],
                    column_config={
                        "duplicates": st.column_config.ListColumn(
                            dups_st[2:],
                            width="large"
                        ),
                    },
                    hide_index=True,
                )
            if non_gene:
                st.data_editor(
                    pd.json_normalize(state.gene_set.validation)['non_hgnc'],
                    column_config={
                        "non_hgnc": st.column_config.ListColumn(
                            non_gene_st[2:],
                            width="large"
                        ),
                    },
                    hide_index=True,
                )
        state.bt_submit_disabled = False


def input_example():
    # Callback because that's the only way it works
    state.gene_set_input = open(f"{ROOT}/code/static/example_gene_list.txt").read()
    state.gene_set_name = "Example gene set"


def main():
    state.bt_submit_disabled = True

    st.sidebar.image(Image.open(f"{ROOT}/code/static/CO_logo_135x72.png"), caption="Code Ocean")
    st.sidebar.title("Enrichment analysis")
    st.sidebar.write(
        """This Streamlit app allows users to submit a list of genes and perform enrichment analysis using Gene Ontology pathways. The app displays enriched pathways and GO terms for the submitted genes, along with relevant statistics such as p-values and FDR corrections."""
    )

    # st.divider()
    st.subheader("Enrichment analysis")
    # col_input, col_genes = st.columns(2)
    # with col_input:
    #     # input_local, input_capsule = st.tabs(["Upload a file", "Browse files in the capsule"])
    #
    #     # with input_local:
    #     uploaded_file = st.file_uploader(
    #         "Upload a set of genes from your local drive",
    #         help="Upload a new line separated set of genes",
    #     )
    #     if uploaded_file is not None:
    #         with open(f"{UPLOAD_FOLDER}/{uploaded_file.name}", "wb") as f:
    #             f.write(uploaded_file.getvalue())
    #         uploaded_value = StringIO(uploaded_file.getvalue().decode("utf-8")).read()
    #         if re.match(hgnc_input_pattern, uploaded_value):
    #             gene_list = uploaded_value
    #             state.bt_submit_disabled = False

    # with col_genes:
    # input_gene_set, settings = st.tabs(["Analyze", "Settings"])
    state.enrich = dict()
    if "results_ready" not in state:
        state.results_ready = False

    analysis, advanced_settings = st.tabs(["Analysis", "Advanced settings"])

    with analysis:
        input_gene_set, settings = st.columns([5, 7])
        with input_gene_set:
            st.text_area(label="Input a set of genes", key="gene_set_input", height=438,
                         placeholder="Input a gene set", label_visibility="collapsed")
            st.text_input(label="Gene set name", key="gene_set_name", placeholder="Input a gene set name",
                          label_visibility="collapsed")

            if "gene_set_input" in state:
                state.bt_submit_disabled = False
                if state.gene_set_input:
                    state.gene_set = GeneSet(state.gene_set_input.split(), state.gene_set_name)

        with settings:
            st.write("Background gene set")
            background_set = st.selectbox("Background gene set",
                                          bg_mapper.keys(),
                                          label_visibility="collapsed")

            st.caption(
                "Specifies the background set of genes. This set validates the input gene set against the chosen organism's genes and serves as a reference for p-value calculations.")
            st.divider()
            st.write('Select libraries')
            libraries = st.multiselect(
                'MSigDB C5 (ontology gene sets)',
                lib_mapper.keys(),
                default=[list(lib_mapper.keys())[0]])

            state.gene_set_libraries = [GeneSetLibrary(f"{ROOT}/data/libraries/{lib_mapper[library]}", name=library) for
                                        library in libraries]
            state.background_gene_set = BackgroundGeneSet(f"{ROOT}/data/backgrounds/{bg_mapper[background_set]}")

        submit, example, placeholder = st.columns([2, 2, 8])
        with submit:
            bt_submit = st.button(
                "Validate and submit", disabled=state.bt_submit_disabled
            )

        with example:
            st.button("Input an example", on_click=input_example)
    with advanced_settings:
        st.write('P-value calculation method')
        st.selectbox('P-value calculation method',
                     options=["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
                     label_visibility="collapsed")
        st.divider()
        st.write("Upload a background gene set")
        bg_custom = st.file_uploader("Upload your background gene set", type=[".txt"], label_visibility="collapsed")
        st.divider()
        st.write("Upload gene set libraries")
        lib_custom = st.file_uploader("Upload gene set libraries", type=[".gmt"], accept_multiple_files=True,
                                      label_visibility="collapsed")

    if bt_submit:
        render_validation()
        # if (state.gene_set_input) and (True in go_libraries.values()):
        if state.gene_set_input:
            n_genes = len(state.gene_set_input.split('\n'))
            if (n_genes <= 10) or (
                    n_genes >= 5000
            ):
                if n_genes <= 100:
                    n_warn = "small"
                elif n_genes >= 2000:
                    n_warn = "big"
                s = 's' if str(n_genes)[-1] != 1 else ''
                st.warning(f"""You've entered {n_genes} gene{s}, which may be {n_warn} and could affect result accuracy. Consider adjusting p-value or log2 Fold Change.  
Estimates for the number of DEGs based on comparison type:
- Similar Conditions (e.g., same cell type, small treatment variations): Dozens to hundreds of DEGs.
- Moderately Different Conditions (e.g., different cell types, moderate drug treatment): Hundreds to thousands.
- Highly Different Conditions (e.g., healthy vs. cancerous tissue): Several thousand DEGs.""")
            with st.spinner("Calculating enrichment"):
                for gene_set_library in state.gene_set_libraries:
                    enrich = Enrichment(state.gene_set, gene_set_library, state.background_gene_set)
                    state.enrich[gene_set_library.name] = enrich
                    with open(f"{ROOT}/results/{enrich.name}.json", "w") as results_snapshot:
                        json.dump(enrich.to_snapshot(), results_snapshot)
                state.results_ready = True
        else:
            if not state.gene_set_input:
                st.error("Please input a newline separated set of genes")
            # if not (True in go_libraries.values()):
            #     st.error("No libraries were selected for the analysis")

    if state.results_ready:
        st.divider()
        for library_name in state.enrich.keys():
            st.subheader(library_name)
            render_results(state.enrich[library_name])

    return


if __name__ == "__main__":
    main()
