import logging
import re
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

import pandas as pd

from src.background_gene_set import BackgroundGeneSet
from src.enrichment import Enrichment
from src.gene_set import GeneSet
from src.gene_set_library import GeneSetLibrary

logger = logging.getLogger(__name__)


class IterativeEnrichment:
    """
    Wrapper for iterative gene set enrichment.

    :param gene_set: List of gene identifiers to analyze.
    :type gene_set: GeneSet
    :param background_gene_set: Background genes file (one gene per line).
    :type background_gene_set: BackgroundGeneSet
    :param gene_set_library: Gene set library GMT file.
    :type gene_set_library: GeneSetLibrary
    :param p_value_method_name: Name of p-value calculation method to pass to Enrichment.
    :type p_value_method_name: str
    :param p_threshold: P-value cutoff for including terms.
    :type p_threshold: float
    :param max_iterations: Maximum number of enrichment iterations to run. None means no limit.
    :type max_iterations: Optional[int]
    """

    def __init__(
        self,
        gene_set: GeneSet,
        gene_set_library: GeneSetLibrary,
        background_gene_set: BackgroundGeneSet,
        p_value_method_name: str = "Fisher's Exact Test",
        name: str = None,
        p_threshold: float = 0.01,
        max_iterations: Optional[int] = None,
    ) -> None:
        self.gene_set = gene_set
        self.gene_set_library = gene_set_library
        self.background_gene_set = background_gene_set
        self.p_value_method_name: str = p_value_method_name
        self.p_threshold: float = p_threshold
        self.max_iterations: Optional[int] = max_iterations
        self.name = (
            name
            if name
            else f"{gene_set.name}_{gene_set_library.name}_{background_gene_set.name}_{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        self._results: List[Dict[str, Any]] = self._compute_enrichment()

    @property
    def results(self) -> List[Dict[str, Any]]:
        """
        The list of iteration records for this enricher.

        :returns: List of dictionaries with keys [iteration, term, library, p-value, genes]
        :rtype: List[Dict[str, Any]]
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
        Perform iterative enrichment, peeling off top terms until no further terms pass p-value threshold.

        :returns: List of iteration records
        :rtype: List[Dict[str, Any]]
        """
        remaining: Set[str] = set(self.gene_set.genes)
        iteration: int = 1
        records: List[Dict[str, Any]] = []

        while True:
            if not remaining:
                logger.info("No genes left; stopping iterative enrichment.")
                break
            if self.max_iterations is not None and iteration > self.max_iterations:
                logger.warning("Reached max_iterations; stopping iterative enrichment.")
                break

            current_set = GeneSet(list(remaining), set(self.background_gene_set.genes))
            try:
                enr = Enrichment(
                    gene_set=current_set,
                    gene_set_library=self.gene_set_library,
                    background_gene_set=self.background_gene_set,
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
                "library": self.gene_set_library.name,
                "p-value": pval,
                "genes": sorted(genes_in_term),
            }
            records.append(record)
            remaining -= genes_in_term
            iteration += 1

        return records

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
        Generate a valid Graphviz DOT for the iterative enrichment network,
        with sanitized, quoted IDs, semicolons, and no duplicates.
        """

        def _sanitize_id(raw: str) -> str:
            """
            Convert raw label into a valid DOT node ID: replace non-alphanumeric with underscores.
            Collapse multiple underscores and strip leading/trailing underscores.
            """
            # replace non-word characters with underscore
            s = re.sub(r"\W+", "_", raw)
            # collapse underscores
            s = re.sub(r"_+", "_", s)
            return s.strip("_")

        nodes: Set[str] = set()
        edges: Set[str] = set()

        # Build nodes and edges
        for rec in self.results:
            term_label = rec.get("term", "")
            # sanitize and quote term ID
            raw_id = f"term_{rec['iteration']}_{term_label}"
            term_id = _sanitize_id(raw_id)
            term_node = (
                f'"{term_id}" '
                f'[label="{term_label} (it {rec["iteration"]})", '
                f'style=filled, fontcolor="white", type="term"];'
            )
            nodes.add(term_node)

            for gene in rec.get("genes", []):
                gene_id = _sanitize_id(f"gene_{gene}")
                gene_node = f'"{gene_id}" [label="{gene}", type="gene"];'
                nodes.add(gene_node)
                # create edge with quoted IDs and semicolon
                edge = f'"{gene_id}" -- "{term_id}";'
                edges.add(edge)

        # Assemble DOT
        lines: List[str] = []
        lines.append("graph iterative_enrichment {")
        lines.append("  graph [layout=neato];")
        lines.append("  node [shape=ellipse];")

        # add nodes
        for node in sorted(nodes):
            lines.append(f"  {node}")
        # add edges
        for edge in sorted(edges):
            lines.append(f"  {edge}")

        lines.append("}")
        return "\n".join(lines)
