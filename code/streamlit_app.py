from math import log10
import json

import plotly.express as px
import streamlit as st
from PIL import Image
import pandas as pd

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from parallel_enrichment import Enrichment
from streamlit import session_state as state

from textarea_onchange import textarea_onchange

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)

ROOT = "/Users/codeocean/PycharmProjects/co_enrichment/"

UPLOAD_FOLDER = f"{ROOT}results/upload"


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


def render_results(result):
    result = result.to_dataframe().head(10)
    result = result.set_index("rank")
    table, bar = st.tabs(["Results", "Bar chart"])
    with table:
        render_table(result)
    with bar:
        render_barchart(result)


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
            st.write(pd.json_normalize(state.gene_set.validation))
        state.bt_submit_disabled = False


def input_example():
    # Callback because that's the only way it works
    state.gene_set_input = open(f"{ROOT}/code/static/example_gene_list.txt").read()
    state.gene_set_name = "Example gene set"


def main():
    state.gene_set_library = GeneSetLibrary(f"{ROOT}/data/libraries/c5.go.cc.v2023.1.Hs.symbols.gmt")
    state.background_gene_set = BackgroundGeneSet(f"{ROOT}/data/backgrounds/hgnc_symbols_2023-01-01.txt")

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
    if "results_ready" not in state:
        state.results_ready = False

    input_gene_set, settings = st.columns([5, 7])
    with input_gene_set:
        textarea_onchange(label="Input a set of genes", key="gene_set_input", height=400, placeholder="Input a gene set", label_visibility="collapsed")
        textarea_onchange(label="Gene set name", key="gene_set_name",height=45, placeholder="Input a gene set name", label_visibility="collapsed")
        submit, example, placholder = st.columns([1, 2, 2])
        if "gene_set_input" in state:
            state.bt_submit_disabled = False
            if state.gene_set_input:
                state.gene_set = GeneSet(state.gene_set_input.split(), state.gene_set_name)
                render_validation()

        with submit:
            bt_submit = st.button(
                "Submit", disabled=state.bt_submit_disabled
            )

        with example:
            st.button("Input an example", on_click=input_example)

        if bt_submit:
            with st.spinner("Calculating enrichment"):
                enrich = Enrichment(state.gene_set, state.gene_set_library, state.background_gene_set)
                with open(f"{ROOT}/results/{enrich.name}.json", "w") as results_snapshot:
                    json.dump(enrich.to_snapshot(), results_snapshot)
                state.results_ready = True

    if state.results_ready:
        st.divider()
        render_results(enrich)

    with settings:
        st.write("Background gene set")
        st.selectbox("Background gene set", options=["HGNC symbols (H. sapiens)", "Mouse"],
                     label_visibility="collapsed")
        st.caption(
            "Specifies the background set of genes. This set validates the input gene set against the chosen organism's genes and serves as a reference for p-value calculations.")
        st.divider()
        st.write('Select libraries')
        libraries = st.multiselect(
            'MSigDB C5 (ontology gene sets)',
            ['GO Biological Process', 'GO Cellular Component', 'GO Molecular Function'],
            default=['GO Biological Process', 'GO Cellular Component', 'GO Molecular Function'])
        st.divider()
        st.write('P-value calculation method')
        st.selectbox('P-value calculation method',
                     options=["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
                     label_visibility="collapsed")

    return


if __name__ == "__main__":
    main()
