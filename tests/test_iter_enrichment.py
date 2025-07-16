from code.background_gene_set import BackgroundGeneSet
from code.gene_set import GeneSet
from code.gene_set_library import GeneSetLibrary
from code.iter_enrichment import IterativeEnrichment
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def create_test_files():
    # Create a background file with genes A, B, C, D
    background_file = ROOT / "tmp" / "background.txt"
    background_file.write_text("A\nB\nC\nD\n")

    # Create a simple GMT library: T1 -> A, B; T2 -> C; T3 -> D
    gmt_file = ROOT / "tmp" / "test_lib.gmt"

    gmt_file.write_text("T1\tdescription\tA\tB\n" "T2\tdescription\tC\n" "T3\tdescription\tD\n")

    return BackgroundGeneSet(str(ROOT / "tmp" / "background.txt")), GeneSetLibrary(str(ROOT / "tmp" / "test_lib.gmt"))


def test_iterative_enrichment_peels_off():
    background_gene_set, gene_set_library = create_test_files()
    genes = GeneSet(["A", "B", "C"])
    enr = IterativeEnrichment(
        gene_set=genes,
        background_gene_set=background_gene_set,
        gene_set_library=gene_set_library,
        p_value_method_name="Fisher's Exact Test",
        p_threshold=1.0,
    )
    df = enr.to_dataframe()

    # Expect two iterations: T1 removes A,B; then T2 removes C
    assert df.shape[0] == 2
    assert df.loc[0, "term"] == "T1"
    assert set(df.loc[0, "genes"]) == {"A", "B"}
    assert df.loc[1, "term"] == "T2"
    assert set(df.loc[1, "genes"]) == {"C"}


def test_to_dot_structure():
    background_gene_set, gene_set_library = create_test_files()
    genes = GeneSet(["A", "B", "C"])
    enr = IterativeEnrichment(
        gene_set=genes,
        background_gene_set=background_gene_set,
        gene_set_library=gene_set_library,
        p_value_method_name="Fisher's Exact Test",
        p_threshold=1.0,
    )
    dot = enr.to_dot()

    # Basic DOT structure checks
    assert dot.startswith("graph iterative_enrichment")
    # Check that term and gene nodes and edges appear
    assert "term_1" in dot
    assert "gene_A -- term_1" in dot
