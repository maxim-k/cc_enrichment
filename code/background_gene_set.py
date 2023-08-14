from typing import List, Set


class BackgroundGeneSet:
    """
    A class to store a set of genes and their size.
    """

    def __init__(self, genes: List[str]) -> None:
        """
        Initialize BackgroundGenes object with a list of genes.

        Args:
            genes: A list of gene names.
        """
        self.genes: Set[str] = set(genes)
        self.size: int = len(self.genes)

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the BackgroundGenes.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
