import logging
from typing import Dict, List, Set

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class GeneSet:
    """
    Class for storing a gene set from a gene list.
    """

    def __init__(
        self,
        gene_list: List[str],
        validation_set: Set,
        name: str = "",
        hgcn: bool = True,
        format: bool = True,
    ) -> None:
        """
        Initialize the class with a list of genes, and two flags - 'hgcn' and 'format'.

        Args:
            gene_list: List of genes
            hgcn: If set to True, removes genes not in 'hgnc_symbols'
            format: If set to True, capitalizes all genes from the input list
        """
        self.genes: Set[str] = set()
        self.name: str = name
        self.size: int = 0
        self.validation: Dict[set[str] : set[str]] = {
            "duplicates": set(),
            "non_valid": set(),
        }
        self.validation_set: Set = validation_set

        if format:
            gene_list = [gene.upper() for gene in gene_list]

        for gene in gene_list:
            if gene in self.genes:
                self.validation["duplicates"].add(gene)
            elif hgcn and gene not in self.validation_set:
                self.validation["non_valid"].add(gene)
            else:
                self.genes.add(gene)

        self.size = len(self.genes)
        logger.info(f"Input Gene Set\n{self.name}\n\t{self.size} genes")
        if self.validation["duplicates"]:
            logger.warning(f"{len(self.validation['duplicates'])} duplicates")
        if self.validation["non_valid"]:
            logger.warning(f"{len(self.validation['non_valid'])} non valid")

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the GeneSet.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
