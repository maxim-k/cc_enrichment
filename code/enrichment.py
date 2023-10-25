import json
from datetime import datetime
import multiprocessing as mp

import pandas as pd

from typing import List, Dict, Any, Tuple
from scipy.stats import fisher_exact
from scipy.stats import hypergeom
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from background_gene_set import BackgroundGeneSet


def compute_pvalue(args: Tuple[GeneSet, BackgroundGeneSet, dict, str]) -> Tuple[str, str, str, List[str], float]:
    """
    Computes the p-value for a given term using Fisher's exact test.
    This function is intended to be used with multiprocessing.Pool.map(),
    which requires functions to take a single argument. Therefore, the
    inputs are passed as a single tuple.

    Args:
        args: A tuple containing the following elements:
            - gene_set (GeneSet): The input gene set
            - background_gene_set (BackgroundGeneSet): The background gene set
            - term (dict): A dictionary representing a term in the gene set library
            - p-value method (str): The name of the method to calculate p-value

    Returns:
        A tuple containing the following elements:
            - term name (str): The name of the term
            - overlap size (int): The size of the overlap between the gene set and the term
            - term description (str): The description of the term
            - overlap genes (list): A list of genes in the overlap
            - p_value (float): The p-value computed by Fisher's exact test
    """
    gene_set, background_gene_set, term, p_value_method_name = args
    term_genes = set(term['genes'])
    n_term_genes = len(term_genes)
    overlap = gene_set.genes & term_genes
    n_overlap = len(overlap)

    # Build contingency table for Fisher's exact test
    contingency_table = [[n_overlap,
                          n_term_genes - n_overlap],
                         [gene_set.size - n_overlap,
                          background_gene_set.size - n_term_genes - gene_set.size + n_overlap]]

    if p_value_method_name == "Fisher's Exact Test":
        _, p_value = fisher_exact(contingency_table)
    elif p_value_method_name == "Chi-squared Test":
        chi2, p_value, _, _ = chi2_contingency(contingency_table)
    elif p_value_method_name == "Hypergeometric Test":
        p_value = hypergeom.sf(n_overlap - 1, background_gene_set.size, n_term_genes, gene_set.size)
    else:
        raise ValueError(f"Unsupported p_value_method: {p_value_method_name}")

    return term['name'], f'{len(overlap)}/{len(term["genes"])}', term['description'], sorted(list(overlap)), p_value


class Enrichment:
    """
    Class for gene set enrichment analysis results.
    """

    def __init__(self, gene_set: GeneSet, gene_set_library: GeneSetLibrary, background_gene_set: BackgroundGeneSet,
                 p_value_method_name="Fisher's Exact Test", name: str = None):
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
        self.p_value_method_name = p_value_method_name
        self.name = name if name else f"{gene_set.name}_{gene_set_library.name}_{background_gene_set.name}_{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
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
        with mp.Pool(mp.cpu_count()) as pool:
            parallel_results = pool.map(compute_pvalue,
                                        [(self.gene_set, self.background_gene_set, term, self.p_value_method_name) for
                                         term in
                                         self.gene_set_library.library])

        # Separate results and p_values for convenience
        p_values = [result[-1] for result in parallel_results]

        # Adjust p-values for multiple testing
        _, p_values_adjusted, _, _ = multipletests(p_values, method='fdr_bh')

        # Rank terms based on their p-values
        ranked_terms = sorted(list(enumerate(parallel_results)), key=lambda x: p_values[x[0]])

        # Format results into a sorted list
        for i, result in ranked_terms:
            term_name, overlap_size, term_description, overlap_genes, _ = result
            results.append({
                'term': term_name,
                'rank': i + 1,
                'description': term_description,
                'overlap': overlap_genes,
                'overlap_size': overlap_size,
                'p-value': p_values[i],
                'fdr': p_values_adjusted[i]
            })

        return results

    def to_dataframe(self):
        """Return the enrichment results as a pandas dataframe."""
        return pd.DataFrame({'rank': [result['rank'] for result in self.results],
                             'term': [result['term'] for result in self.results],
                             'overlap': [result['overlap'] for result in self.results],
                             'overlap_size': [result['overlap_size'] for result in self.results],
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

    def to_snapshot(self) -> Dict:
        """Return the snapshot of input parameters and the enrichment results as a JSON string."""
        return {
            "input_gene_set": list(self.gene_set.genes),
            "background": self.background_gene_set.name,
            self.gene_set_library.name: self.results
        }
