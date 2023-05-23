import json

import pandas as pd

from typing import List, Dict, Any
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from background_gene_set import BackgroundGeneSet


class Enrichment:
    """
    Class for gene set enrichment analysis results.
    """

    def __init__(self, gene_set: GeneSet, gene_set_library: GeneSetLibrary, background_gene_set: BackgroundGeneSet):
        """
        Initialize the class with gene set, gene set library, and background gene set.

        Args:
            gene_set: Input gene set
            gene_set_library: Gene set library
            background_gene_set: Background gene set
        """
        self.gene_set = gene_set
        self.gene_set_library = gene_set_library
        self.background_gene_set = background_gene_set
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        Getter for _results.

        Returns:
            A list containing dictionaries of enrichment results
        """
        return self._results

    @results.setter
    def results(self, value: List[Dict[str, Any]]) -> None:
        """
        Setter for _results.

        Args:
            value: A list containing dictionaries of enrichment results
        """
        self._results = value

    def _compute_enrichment(self) -> List[Dict[str, Any]]:
        """
        Computes gene set enrichment analysis.

        Returns:
            A list containing dictionaries of enrichment results
        """
        results = []
        p_values = []
        for term in self.gene_set_library.library:
            term_genes = set(term['genes'])
            overlap = self.gene_set.genes & term_genes

            # Build contingency table for Fisher's exact test
            contingency_table = [[len(overlap),
                                  len(term_genes) - len(overlap)],
                                 [len(self.gene_set.genes) - len(overlap),
                                  len(self.background_gene_set.genes) - len(term_genes) - len(self.gene_set.genes) + len(overlap)]]

            # Perform Fisher's exact test
            _, p_value = fisher_exact(contingency_table)
            p_values.append(p_value)

        # Adjust p-values for multiple testing
        _, p_values_adjusted, _, _ = multipletests(p_values, method='fdr_bh')

        # Rank terms based on their p-values
        ranked_terms = sorted(list(enumerate(self.gene_set_library.library)), key=lambda x: p_values[x[0]])

        for i, term in ranked_terms:
            results.append({
                'term': term['name'],
                'rank': i + 1,
                'description': term['description'],
                'overlap': sorted(list(self.gene_set.genes & set(term['genes']))),
                'p-value': p_values[i],
                'fdr': p_values_adjusted[i]
            })
        return results

    def to_dataframe(self):
        """Return the enrichment results as a pandas dataframe."""
        return pd.DataFrame({'rank': [result['rank'] for result in self.results],
                             'term': [result['term'] for result in self.results],
                             'overlap': [result['overlap'] for result in self.results],
                             'p-value': [result['p-value'] for result in self.results],
                             'fdr': [result['fdr'] for result in self.results]
                             })

    def to_json(self):
        """Return the enrichment results as a JSON string."""
        return json.dumps(self.results, indent=4, separators=(',', ': '))

    def to_html(self):
        """Return the enrichment results as an HTML page."""
        return self.to_dataframe().to_html()

    def to_tsv(self):
        """Return the enrichment results as a TSV spreadsheet."""
        return self.to_dataframe().to_csv(sep='\t')
