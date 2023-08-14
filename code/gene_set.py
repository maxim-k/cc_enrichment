from typing import List, Set, Dict

ROOT = "/Users/codeocean/PycharmProjects/co_enrichment/"
class GeneSet:
    """
    Class for storing a gene set from a gene list.
    """
    with open(f'{ROOT}/data/hgnc_symbols_2023-01-01.txt', 'r') as f:
        hgnc_symbols: Set[str] = set(line.strip() for line in f)

    def __init__(self, gene_list: List[str], hgcn: bool = True, format: bool = True) -> None:
        """
        Initialize the class with a list of genes, and two flags - 'hgcn' and 'format'.

        Args:
            gene_list: List of genes
            hgcn: If set to True, removes genes not in 'hgnc_symbols'
            format: If set to True, capitalizes all genes from the input list
        """
        self.genes: Set[str] = set()
        self.size: int = 0
        self.validation: Dict[set[str]: set[str]] = {'duplicates': set(), 'non_hgnc': set()}

        if format:
            gene_list = [gene.upper() for gene in gene_list]

        for gene in gene_list:
            if gene in self.genes:
                self.validation['duplicates'].add(gene)
            elif hgcn and gene not in self.hgnc_symbols:
                self.validation['non_hgnc'].add(gene)
            else:
                self.genes.add(gene)

        self.size = len(self.genes)

    def has_gene(self, gene: str) -> bool:
        """
        Check if the given gene is present in the GeneSet.

        Args:
            gene: A gene name.

        Returns:
            True if the gene is present, False otherwise.
        """
        return gene in self.genes
