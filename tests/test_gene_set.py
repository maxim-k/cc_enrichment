from code.gene_set import GeneSet


class TestGeneSet:

    def test_gene_set_with_defaults(self):
        gs = GeneSet(['A1bg', 'bgn', 'BGN', 'celf2-as1', 'ddx50p1', 'non_hgnc', 'A1bg'])
        assert gs.size == 4
        assert gs.validation == {'duplicates': {'A1BG', 'BGN'}, 'non_hgnc': {'NON_HGNC'}}

    def test_gene_set_without_hgcn(self):
        gs = GeneSet(['A1bg', 'bgn', 'BGN', 'celf2-as1', 'ddx50p1', 'non_hgnc', 'A1bg'], hgcn=False)
        assert gs.size == 5
        assert gs.validation == {'duplicates': {'A1BG', 'BGN'}, 'non_hgnc': set()}

    def test_gene_set_without_format(self):
        gs = GeneSet(['A1bg', 'bgn', 'BGN', 'celf2-as1', 'ddx50p1', 'non_hgnc', 'A1bg'], format=False)
        assert gs.size == 1
        assert gs.validation == {'duplicates': set(), 'non_hgnc': {'A1bg', 'bgn', 'celf2-as1', 'ddx50p1', 'non_hgnc'}}
