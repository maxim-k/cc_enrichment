import re
from io import StringIO
from math import log10

import plotly.express as px
import streamlit as st
from PIL import Image

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from parallel_enrichment import Enrichment

# from treeview_streamlit_component import treeview_streamlit_component

st.set_page_config(
    page_title="Enrichment Analysis", layout="wide", initial_sidebar_state="expanded"
)

UPLOAD_FOLDER = "results/upload"


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


def input_example():
    # Callback because that's the only way it works
    st.session_state.text_input = open("static/example_gene_list.txt").read()


def main():
    gene_list = ""
    gene_set_library = GeneSetLibrary("../data/c5.go.bp.v2023.1.Hs.symbols.gmt")
    background_gene_set = BackgroundGeneSet(
        open("../data/hgnc_symbols_2023-01-01.txt", "r").read().split()
    )
    hgnc_input_pattern = r"^(?:\b\w{1,30}\b\n?)+$"
    st.session_state.bt_submit_disabled = True

    st.sidebar.image(Image.open("static/CO_logo_135x72.png"), caption="Code Ocean")
    st.sidebar.title("Enrichment analysis")
    st.sidebar.write(
        """This Streamlit app allows users to submit a list of genes and perform enrichment analysis using Gene Ontology pathways. The app displays enriched pathways and GO terms for the submitted genes, along with relevant statistics such as p-values and FDR corrections."""
    )

    st.divider()
    st.subheader("Input a gene set")
    col_input, col_genes = st.columns(2)
    with col_input:
        # input_local, input_capsule = st.tabs(["Upload a file", "Browse files in the capsule"])

        # with input_local:
        uploaded_file = st.file_uploader(
            "Upload a set of genes from your local drive",
            help="Upload a new line separated set of genes",
        )
        if uploaded_file is not None:
            with open(f"{UPLOAD_FOLDER}/{uploaded_file.name}", "wb") as f:
                f.write(uploaded_file.getvalue())
            uploaded_value = StringIO(uploaded_file.getvalue().decode("utf-8")).read()
            if re.match(hgnc_input_pattern, uploaded_value):
                gene_list = uploaded_value
                st.session_state.bt_submit_disabled = False

    with col_genes:
        st.text_area(
            label="Input a set of genes", key="text_input", value=gene_list, height=400
        )
        if re.match(hgnc_input_pattern, st.session_state.text_input):
            gene_list = st.session_state.text_input
            st.session_state.bt_submit_disabled = False

        submit, example, blank = st.columns([3, 4, 3])
        with submit:
            bt_submit = st.button(
                "Submit", disabled=st.session_state.bt_submit_disabled
            )

        with example:
            st.button("Input an example", on_click=input_example)

    if st.session_state.text_input:
        if bt_submit:
            gene_set = GeneSet(st.session_state.text_input.split())
            with st.spinner("Calculating enrichment"):
                enrich = Enrichment(gene_set, gene_set_library, background_gene_set)

            render_results(enrich)

    return


if __name__ == "__main__":
    main()
