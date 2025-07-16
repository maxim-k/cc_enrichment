import logging
from typing import List, Dict, Any, Set, Optional

import pandas as pd

from background_gene_set import BackgroundGeneSet
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from enrichment import Enrichment

logger = logging.getLogger(__name__)


class IterativeEnricher:
    """
    Wrapper for iterative gene set enrichment.

    :param gene_list: List of gene identifiers to analyze.
    :type gene_list: List[str]
    :param background_file: Path to the background genes file (one gene per line).
    :type background_file: str
    :param gmt_file: Path to the gene set library GMT file.
    :type gmt_file: str
    :param p_value_method_name: Name of p-value calculation method to pass to Enrichment.
    :type p_value_method_name: str
    :param p_threshold: P-value cutoff for including terms.
    :type p_threshold: float
    :param max_iterations: Maximum number of enrichment iterations to run. None means no limit.
    :type max_iterations: Optional[int]
    """

    def __init__(
        self,
        gene_list: List[str],
        background_file: str,
        gmt_file: str,
        p_value_method_name: str = "Fisher's Exact Test",
        p_threshold: float = 0.01,
        max_iterations: Optional[int] = None,
    ) -> None:
        self.background: BackgroundGeneSet = BackgroundGeneSet(background_file)
        self.library: GeneSetLibrary = GeneSetLibrary(gmt_file)
        self.master_genes: GeneSet = GeneSet(gene_list, set(self.background.genes))
        self.p_value_method_name: str = p_value_method_name
        self.p_threshold: float = p_threshold
        self.max_iterations: Optional[int] = max_iterations
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    def _compute_enrichment(self) -> List[Dict[str, Any]]:
        """
        Perform iterative enrichment, peeling off top terms until no further terms pass p-value threshold.

        :returns: List of iteration records
        :rtype: List[Dict[str, Any]]
        """
        remaining: Set[str] = set(self.master_genes.genes)
        iteration: int = 1
        records: List[Dict[str, Any]] = []

        while True:
            if not remaining:
                logger.info("No genes left; stopping iterative enrichment.")
                break
            if self.max_iterations is not None and iteration > self.max_iterations:
                logger.warning("Reached max_iterations; stopping iterative enrichment.")
                break

            current_set = GeneSet(list(remaining), set(self.background.genes))
            try:
                enr = Enrichment(
                    gene_set=current_set,
                    gene_set_library=self.library,
                    background_gene_set=self.background,
                    p_value_method_name=self.p_value_method_name,
                )
            except Exception as e:
                logger.error(f"Enrichment failed at iteration {iteration}: {e}")
                break

            results = enr.results
            if not results:
                logger.info("No enrichment results; terminating.")
                break

            top = results[0]
            pval = top.get("p-value")
            if pval is None or pval >= self.p_threshold:
                logger.info("Top term p-value >= threshold; terminating.")
                break

            genes_in_term: Set[str] = set(top.get("overlap", []))
            record: Dict[str, Any] = {
                "iteration": iteration,
                "term": top.get("term", ""),
                "library": self.library.name,
                "p-value": pval,
                "genes": sorted(genes_in_term),
            }
            records.append(record)
            remaining -= genes_in_term
            iteration += 1

        return records

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        The list of iteration records for this enricher.

        :returns: List of dictionaries with keys [iteration, term, library, p-value, genes]
        :rtype: List[Dict[str, Any]]
        """
        return self._results

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert iteration records to a pandas DataFrame.

        :returns: DataFrame of iteration records
        :rtype: pandas.DataFrame
        """
        return pd.DataFrame(self.results)

    def to_tsv(self) -> str:
        """
        Serialize iteration records as a TSV string.

        :returns: TSV-formatted string
        :rtype: str
        """
        return self.to_dataframe().to_csv(sep="\t", index=False)

    def to_json(self) -> str:
        """
        Serialize iteration records as a JSON-formatted string.

        :returns: JSON-formatted string
        :rtype: str
        """
        import json

        return json.dumps(self.results, indent=2)

    def to_dot(self) -> str:
        """
        Generate a Graphviz DOT representation of the iterative enrichment network.

        :returns: DOT-format string
        :rtype: str
        """
        # Define a simple color palette
        palette: List[str] = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
        lib_color = palette[0]

        term_nodes: Set[str] = set()
        gene_nodes: Set[str] = set()
        edges: List[str] = []

        for rec in self.results:
            term_id = f"term_{rec['iteration']}"
            term_label = rec.get("term", "")
            term_nodes.add(f'{term_id} [label="{term_label}", color="{lib_color}"]')
            for gene in rec.get("genes", []):
                gene_id = f"gene_{gene}"
                gene_nodes.add(f'{gene_id} [label="{gene}"]')
                edges.append(f"{gene_id} -- {term_id}")

        lines: List[str] = [
            "graph iterative_enrichment {",
            "  graph [layout=neato];",
        ]
        for node in term_nodes | gene_nodes:
            lines.append(f"  {node}")
        for edge in edges:
            lines.append(f"  {edge}")
        lines.append("}")

        return "\n".join(lines)
