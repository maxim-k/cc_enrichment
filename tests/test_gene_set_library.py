from gene_set_library import GeneSetLibrary


class TestGeneSetLibrary:
    def setup_method(self):
        self.gmt_file_path = 'sample.gmt'
        self.library = GeneSetLibrary(self.gmt_file_path)

    def test_num_terms(self):
        assert self.library.num_terms == 5

    def test_unique_genes(self):
        assert len(self.library.unique_genes) == 171
