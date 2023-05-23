from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from timeit import default_timer as timer


class TestEnrichment:

    def test_enrichment_results(self):
        gene_set = GeneSet(
            ['ACO1', 'ADH5', 'ADK', 'BRI3', 'C1D', 'CAT', 'CNO', 'DDT', 'ELL3', 'ENY2', 'ESM1', 'FAH', 'GBE1', 'GK5',
             'GLO1', 'GYK', 'HPN', 'HYI', 'IPP', 'KLF1', 'LIFR', 'LYRM2', 'LYRM5', 'MCAT', 'MPP7', 'MUT', 'MYNN',
             'MYO6', 'NEO1', 'NOL7', 'NPY', 'NSUN3', 'NUPL2', 'OS'])
        gene_set_library = GeneSetLibrary('./tests/sample.gmt')
        background_gene_set = BackgroundGeneSet(open('./hgnc_symbols_2023-01-01.txt', 'r').read().split())
        start = timer()
        enrich = Enrichment(gene_set, gene_set_library, background_gene_set)

        results = enrich.results
        end = timer()
        print(end - start)
        print(enrich.to_json())
        assert isinstance(results, list)
        for result in results:
            assert isinstance(result, dict)
            assert isinstance(result['term'], str)
            assert isinstance(result['rank'], int)
            assert isinstance(result['description'], str)
            assert isinstance(result['overlap'], list)
            assert isinstance(result['p-value'], float)
            assert isinstance(result['fdr'], float)
