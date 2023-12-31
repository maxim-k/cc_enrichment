from pathlib import Path
from code.background_gene_set import BackgroundGeneSet

ROOT = Path(__file__).resolve().parent.parent


class TestBackgroundGeneSet:
    def test_background_genes_creation(self):
        genes_list = "A1BG\nBGN\nCELF2-AS1"
        bg = BackgroundGeneSet(genes_list)
        assert bg.size == 3
        assert bg.has_gene("A1BG")
        assert bg.has_gene("BGN")
        assert bg.has_gene("CELF2-AS1")

    def test_background_genes_empty_creation(self):
        bg = BackgroundGeneSet("")
        assert bg.size == 0
        assert not bg.has_gene("A1BG")
        assert not bg.has_gene("BGN")
        assert not bg.has_gene("CELF2-AS1")

    def test_background_genes_with_duplicates(self):
        genes_list = "A1BG\nBGN\nA1BG"
        bg = BackgroundGeneSet(genes_list)
        assert bg.size == 2
        assert bg.has_gene("A1BG")
        assert bg.has_gene("BGN")
