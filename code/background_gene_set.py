from typing import Set
from pathlib import Path


class BackgroundGeneSet:
    """
    A class to store a set of genes and their size.
    """

    def __init__(self, background_file_path: str, name: str = "", organism: str = "Homo Sapiens") -> None:
        """
        Initialize BackgroundGenes object with a list of genes.

        Args:
            genes: A list of gene names.
        """
        self.genes: Set[str] = self._load_from_file(background_file_path)
        self.size: int = len(self.genes)
        self.name = name if name else Path(background_file_path).stem
        self.organism = organism

    def _load_from_file(self, background_file_path: str) -> Set[str]:
        """
        Load library from a .gmt file
        Args:
            background_file_path: Path to the background file

        Returns:
            Set of genes representing the background
        """
        return set(open(background_file_path, "r").read().split('\n'))

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the BackgroundGenes.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
