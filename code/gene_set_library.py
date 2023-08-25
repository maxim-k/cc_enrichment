from typing import List, Dict, Set
from pathlib import Path


class GeneSetLibrary:
    """
    Class to represent gene set library.
    """

    def __init__(self, gmt_file_path: str, name: str = "", organism: str = "Homo Sapiens") -> None:
        """
        Initialize gene set library from a .gmt file

        :param gmt_file_path: Path to the .gmt file
        """
        self.library = self._load_from_gmt(gmt_file_path)
        self.num_terms = len(self.library)
        self.unique_genes = self._compute_unique_genes()
        self.size = len(self.unique_genes)
        self.name = name if name else Path(gmt_file_path).stem
        self.organism = organism

    def _load_from_gmt(self, gmt_file_path: str) -> List[Dict[str, List[str]]]:
        """
        Load library from a .gmt file
        Args:
            gmt_file_path: Path to the .gmt file

        Returns:
            List of dictionaries representing the library
        """
        library = []
        with open(gmt_file_path, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')
                term = {
                    'name': parts[0],
                    'description': parts[1],
                    'genes': parts[2:]
                }
                library.append(term)
        return library

    def _compute_unique_genes(self) -> Set[str]:
        """
        Compute the set of unique genes in the library.

        :return: Set of unique genes
        """
        unique_genes = set()
        for term in self.library:
            unique_genes.update(term['genes'])
        return unique_genes

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the GeneSet.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.unique_genes
