import json
import logging
from io import StringIO
from pathlib import Path
from typing import Dict, List

import math
import streamlit as st
from PIL import Image
from streamlit import session_state as state

# Existing local imports (kept exactly as in original file)
from background_gene_set import BackgroundGeneSet
from src.enrichment import Enrichment
from src.gene_set import GeneSet
from src.gene_set_library import GeneSetLibrary
from src.iter_enrichment import IterativeEnrichment  # NEW
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


# ---------- Helper functions for NEW iterative mode ----------

def _ensure_base_state():
    """Initialize shared state keys once (idempotent)."""
    if "enrich" not in state:
        state.enrich = dict()
    if "iter_results" not in state:
        # { library_name: List[iteration dicts] }
        state.iter_results: Dict[str, List[dict]] = {}
    if "iter_graph_parts" not in state:
        # store per-library DOT snippets (term + gene nodes + edges)
        state.iter_graph_parts: Dict[str, dict] = {}
    if "results_ready" not in state:
        state.results_ready = False
    if "iter_ready" not in state:
        state.iter_ready = False


def _build_iterative_tables_download(all_iter_results: Dict[str, List[dict]]) -> str:
    """
    Combine per-library iterative records into a single TSV string.
    Columns: library, iteration, term, p-value, genes (comma joined), -log10_p.
    """
    rows: List[str] = ["library\titeration\tterm\tp-value\t-log10_p\tgenes"]
    for lib, records in all_iter_results.items():
        for rec in records:
            p = rec.get("p-value", float("nan"))
            rows.append(
                f"{lib}\t{rec.get('iteration')}\t{rec.get('term','')}\t{p}"
                f"\t{(-math.log10(p) if (p and p>0) else '')}\t{','.join(rec.get('genes', []))}"
            )
    return "\n".join(rows)


def _merge_iterative_dot(per_lib_parts: Dict[str, dict]) -> str:
    """
    Merge per-library node/edge parts into one DOT graph with legend.
    Each entry in per_lib_parts: { 'nodes_term': set(str), 'nodes_gene': set(str), 'edges': set(str), 'color': str }.
    """
    lines: List[str] = [
        "graph iterative_enrichment_all {",
        "  graph [layout=neato, overlap=false];",
        "  node [shape=ellipse, fontname=\"Helvetica\"];",
    ]
    # Legend cluster
    lines.append("  subgraph cluster_legend {")
    lines.append("    label=\"Legend\";")
    legend_y = 0
    for i, (lib, parts) in enumerate(per_lib_parts.items(), start=1):
        color = parts["color"]
        # tiny invisible nodes for legend
        lines.append(
            f'    legend_{i} [label="{lib}", shape=box, style=filled, fillcolor="{color}", fontcolor="white"];'
        )
        legend_y += 1
    lines.append("  }")

    # Add nodes / edges
    for parts in per_lib_parts.values():
        for n in parts["nodes_term"]:
            lines.append(f"  {n}")
        for n in parts["nodes_gene"]:
            lines.append(f"  {n}")
        for e in parts["edges"]:
            lines.append(f"  {e}")

    lines.append("}")
    return "\n".join(lines)


def _extract_graph_parts(library_name: str, records: List[dict], color: str) -> dict:
    """
    Build node/edge fragments for later merged DOT. Distinguish term nodes by iteration+library.
    Gene node shared across libraries -> consistent id gene_<gene>.
    """
    nodes_term = set()
    nodes_gene = set()
    edges = set()
    for rec in records:
        term_id = f"term_{library_name}_{rec['iteration']}"
        term_label = rec.get("term", "")
        nodes_term.add(
            f'{term_id} [label="{term_label}\\n(it {rec["iteration"]})", '
            f'style=filled, fillcolor="{color}", fontcolor="white"]'
        )
        for gene in rec.get("genes", []):
            gene_id = f"gene_{gene}"
            nodes_gene.add(f'{gene_id} [label="{gene}"]')
            edges.add(f"{gene_id} -- {term_id}")
    return {
        "nodes_term": nodes_term,
        "nodes_gene": nodes_gene,
        "edges": edges,
        "color": color,
    }


# Static palette (reuse from iterative module base but extended)
_PALETTE = [
    "#1b9e77",
    "#d95f02",
    "#7570b3",
    "#e7298a",
    "#66a61e",
    "#e6ab02",
    "#a6761d",
    "#666666",
]


# ---------- Main App ----------

def main() -> None:
    """
    Streamlit application entry point.
    Regular and Iterative Enrichment modes supported.
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

    _ensure_base_state()

    # ---------- Mode Switch ----------
    mode = st.radio(
        "Mode",
        options=["Regular", "Iterative"],
        horizontal=True,
        help="Choose 'Regular' for standard enrichment with FDR; 'Iterative' to peel top terms iteratively.",
        key="analysis_mode",
    )

    st.subheader(f"Enrichment analysis — {mode} mode")

    # Shared consistent initial mappings each run
    state.lib_mapper = update_aliases("libraries")
    state.bg_mapper = update_aliases("backgrounds")
    state.advanced_settings_changed = False
    state.bt_submit_disabled = True

    analysis, advanced_settings = st.tabs(["Analysis", "Advanced settings"])

    # ---------- ANALYSIS TAB (common inputs) ----------
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

            st.selectbox(
                "Or select a file from the `data` folder",
                ["Select ..."] + gene_set_files,
                index=0,
                on_change=update_text_widgets,
                key="selected_file",
            )

    with settings:
        # Background
        state.background_set = st.selectbox(
            "Background gene set", state.bg_mapper.keys()
        )
        st.caption(
            "Specifies the background set of genes. This set validates the input gene set against the chosen organism's genes and serves as a reference for p-value calculations."
        )

        # Library selection
        state.libraries = st.multiselect(
            "Select libraries", state.lib_mapper.keys(), default=None
        )

        # Construct library objects
        if ("libraries" in state) and ("lib_mapper" in state) and state.libraries:
            state.gene_set_libraries = [
                GeneSetLibrary(
                    str(ROOT / "data" / "libraries" / state.lib_mapper[library]),
                    name=library,
                )
                for library in state.libraries
            ]
        else:
            state.gene_set_libraries = []

        # Background object + GeneSet
        if ("background_set" in state) and ("bg_mapper" in state):
            state.background_gene_set = BackgroundGeneSet(
                str(
                    ROOT
                    / "data"
                    / "backgrounds"
                    / state.bg_mapper[state.background_set]
                )
            )
            if "gene_set_input" in state and state.gene_set_input:
                state.gene_set = GeneSet(
                    state.gene_set_input.split(),
                    state.background_gene_set.genes,
                    state.gene_set_name,
                )

        # Iterative-specific controls (threshold etc.)
        if mode == "Iterative":
            st.markdown("**Iterative parameters**")
            state.iter_p_threshold = st.number_input(
                "P-value threshold",
                min_value=1e-10,
                max_value=0.5,
                value=0.01,
                step=0.001,
                format="%.4f",
                help="Stop after top term p-value >= threshold.",
            )
            state.iter_max_iter = st.number_input(
                "Max iterations (0 = no limit)",
                min_value=0,
                max_value=500,
                value=0,
                step=1,
                help="Optional safety cap on iterations.",
            )

    # ---------- ACTION BUTTONS ----------
    submit_col, example_col, _ = st.columns([9, 8, 29])

    # Common enable condition (gene set + background + libraries)
    ready_common = all(
        key in state and state[key]
        for key in ["gene_set", "background_gene_set", "gene_set_libraries"]
    )

    with submit_col:
        if mode == "Regular":
            state.bt_submit_disabled = not ready_common
            bt_submit = st.button(
                "Validate and submit", disabled=state.bt_submit_disabled, key="bt_regular"
            )
        else:
            state.bt_iter_disabled = not ready_common
            bt_iter_submit = st.button(
                "Run iterative enrichment",
                disabled=state.bt_iter_disabled,
                key="bt_iterative",
            )

    with example_col:
        st.button("Input an example", on_click=input_example)

    # ---------- ADVANCED SETTINGS TAB ----------
    with advanced_settings:
        if mode == "Regular":
            n_results = st.slider(
                "Number of results to display",
                min_value=1,
                max_value=100,
                value=10,
                step=1,
                key="n_results_regular",
            )
        else:
            # still define but hidden for iterative result reuse (we could allow controlling display depth)
            n_results = st.slider(
                "Number of results to display (regular mode only)",
                min_value=1,
                max_value=100,
                value=10,
                step=1,
                disabled=True,
                key="n_results_iter_hidden",
            )

        state.p_val_method = st.selectbox(
            "P-value calculation method",
            options=["Fisher's Exact Test", "Hypergeometric Test", "Chi-squared Test"],
            key="p_val_method_select",
        )
        if state.p_val_method != "Fisher's Exact Test":
            state.advanced_settings_changed = True

        # Background upload
        state.bg_custom = st.file_uploader(
            "Upload your background gene set", type=[".txt"]
        )
        if state.bg_custom is not None:
            bg_file = (ROOT / "data" / "backgrounds" / state.bg_custom.name).open("wb")
            bg_file.write(state.bg_custom.getvalue())
            state.advanced_settings_changed = True

        # Library upload(s)
        state.libs_custom = st.file_uploader(
            "Upload gene set libraries",
            type=[".gmt"],
            accept_multiple_files=True,
            on_change=update_aliases,
            args=("libraries",),
        )
        if state.libs_custom:
            for lib_custom in state.libs_custom:
                lib_file = (ROOT / "data" / "libraries" / lib_custom.name).open("wb")
                lib_file.write(lib_custom.getvalue())
                state.advanced_settings_changed = True

        # Apply settings button
        if state.advanced_settings_changed:
            if st.button("Apply settings"):
                logger.info("Running with custom settings")
                st.success("Settings applied")
        else:
            with st.empty():
                st.button("Apply settings", disabled=True)

    # ---------- REGULAR MODE EXECUTION ----------
    if mode == "Regular" and 'bt_submit' in locals() and bt_submit:
        logger.info("Validating and submitting genes for enrichment analysis (regular)")
        render_validation()
        if state.gene_set_input and ready_common:
            n_genes = len(state.gene_set_input.split("\n"))
            if (n_genes <= 100) or (n_genes >= 5000):
                n_warn = "small" if n_genes <= 100 else "big"
                s = "s" if str(n_genes)[-1] != "1" else ""
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
                logger.info(
                    f"Enrichment results for {gene_set_library.name} are ready (regular)"
                )
                state.results_ready = True
        else:
            if not state.gene_set_input:
                logger.error("Please input a newline separated set of genes")
                st.error("Please input a newline separated set of genes")
            if not state.gene_set_libraries:
                logger.error("No libraries were selected for the analysis")
                st.error("No libraries were selected for the analysis")
            if not getattr(state, "background_gene_set", None):
                logger.error("No background gene set was selected for the analysis")
                st.error("No background gene set was selected for the analysis")

    if mode == "Regular" and state.results_ready:
        logger.info("Displaying enrichment results (regular)")
        st.markdown(
            f'Download all results as {download_link(collect_results(state.enrich), "results", "tsv")}',
            unsafe_allow_html=True,
        )
        for library_name in state.enrich.keys():
            render_results(state.enrich[library_name], library_name, n_results)
        state.results_ready = False

    # ---------- ITERATIVE MODE EXECUTION ----------
    if mode == "Iterative" and 'bt_iter_submit' in locals() and bt_iter_submit:
        logger.info("Running iterative enrichment")
        # Validation similar to regular
        if not state.gene_set_input:
            st.error("Please input a newline separated set of genes")
        if not state.gene_set_libraries:
            st.error("No libraries were selected for the analysis")
        if not getattr(state, "background_gene_set", None):
            st.error("No background gene set was selected for the analysis")

        if ready_common and state.gene_set_input:
            state.iter_results.clear()
            state.iter_graph_parts.clear()
            with st.spinner("Running iterative enrichment"):
                for idx, gene_set_library in enumerate(state.gene_set_libraries):
                    color = _PALETTE[idx % len(_PALETTE)]
                    logger.info(
                        f"Iterative enrichment for {gene_set_library.name} (color {color})"
                    )
                    it_enr = IterativeEnrichment(
                        gene_set=state.gene_set,
                        gene_set_library=gene_set_library,
                        background_gene_set=state.background_gene_set,
                        p_value_method_name=state.p_val_method,
                        p_threshold=state.iter_p_threshold,
                        max_iterations=(
                            None if state.iter_max_iter == 0 else state.iter_max_iter
                        ),
                    )
                    recs = it_enr.results
                    state.iter_results[gene_set_library.name] = recs
                    # Build partial graph parts for merged DOT
                    state.iter_graph_parts[gene_set_library.name] = _extract_graph_parts(
                        gene_set_library.name, recs, color
                    )
            state.iter_ready = True

    # ---------- ITERATIVE MODE RESULT RENDERING ----------
    if mode == "Iterative" and state.iter_ready:
        # Sub-tabs: Iterations, Network
        iter_tab, network_tab = st.tabs(["Iterations", "Network"])

        # Download (combined TSV)
        combined_tsv = _build_iterative_tables_download(state.iter_results)
        st.markdown(
            f'Download iterative results as {download_link(combined_tsv, "iterative_results", "tsv")}',
            unsafe_allow_html=True,
        )

        with iter_tab:
            if not state.iter_results:
                st.info("No iterative enrichment results.")
            else:
                for lib_name, records in state.iter_results.items():
                    st.markdown(f"### {lib_name}")
                    if not records:
                        st.write("_No terms passed threshold._")
                        continue
                    # Build dataframe manually (avoid external dependency)
                    import pandas as pd  # local import to mirror pattern
                    df_rows = []
                    for rec in records:
                        p = rec.get("p-value")
                        df_rows.append(
                            {
                                "Iteration": rec["iteration"],
                                "Term": rec["term"],
                                "p-value": p,
                                "-log10(p)": (-math.log10(p) if p and p > 0 else None),
                                "Genes removed": ", ".join(rec.get("genes", [])),
                            }
                        )
                    df = pd.DataFrame(df_rows)
                    st.dataframe(df, use_container_width=True)

                    # Bar chart of –log10(p) vs Iteration
                    # Using basic altair via st.bar_chart for simplicity
                    chart_df = df[["Iteration", "-log10(p)"]].set_index("Iteration")
                    st.bar_chart(chart_df)

                    # Per-library TSV download
                    per_lib_tsv_io = StringIO()
                    per_lib_tsv_io.write(
                        "iteration\tterm\tp-value\tgenes\n"
                    )
                    for rec in records:
                        per_lib_tsv_io.write(
                            f"{rec['iteration']}\t{rec['term']}\t{rec['p-value']}\t{','.join(rec['genes'])}\n"
                        )
                    per_lib_tsv = per_lib_tsv_io.getvalue()
                    st.markdown(
                        f'Download {lib_name} results as {download_link(per_lib_tsv, f"iter_{lib_name}", "tsv")}',
                        unsafe_allow_html=True,
                    )

        with network_tab:
            if not state.iter_results:
                st.info("No network to display.")
            else:
                merged_dot = _merge_iterative_dot(state.iter_graph_parts)
                st.graphviz_chart(merged_dot)
                st.markdown(
                    f'Download merged network DOT as {download_link(merged_dot, "iterative_network", "dot")}',
                    unsafe_allow_html=True,
                )
        # Reset flag so repeated reruns don't duplicate unless button pressed
        state.iter_ready = False

    logger.info("Finishing the Streamlit app")
    return


if __name__ == "__main__":
    main()
