import typer
from typing import List, Optional

from background_gene_set import BackgroundGeneSet
from enrichment import Enrichment
from gene_set import GeneSet
from gene_set_library import GeneSetLibrary
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
app = typer.Typer()


def run_enrichment():
    return
@app.command()
def main(
    gene_sets: List[Path] = typer.Option(
        None,
        "--gene-sets", "-g",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Paths to gene set files."
    ),
    background: Optional[Path] = typer.Option(
        None,
        "--background", "-b",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Path to the background gene set file."
    ),
    libraries: List[Path] = typer.Option(
        None,
        "--libraries", "-l",
        exists=True,
        file_okay=True,
        dir_okay=False,
        help="Paths to gene set library files."
    ),

    p_value_method: str = typer.Option(
        "fishers_exact",
        "--method", "-m",
        help="P-value calculation method."
    )
):
    # Default values handling
    if not gene_sets:
        gene_sets = list((ROOT / "data/gene_lists").glob("*.txt"))

    if not background:
        background_files = list((ROOT / "data/backgrounds").glob("*.txt"))
        background = background_files[0] if background_files else None

    if not libraries:
        libraries = list((ROOT / "data/libraries").glob("*.gmt"))

    # Ensure that the resulting variables are not empty
    if not gene_sets or not libraries:
        typer.echo("Error: Gene sets and libraries cannot be empty.", err=True)
        raise typer.Exit(code=1)

    run_enrichment()


    return


if __name__ == "__main__":
    main()
