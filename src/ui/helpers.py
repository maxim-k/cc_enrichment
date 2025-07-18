import logging
from pathlib import Path

from streamlit import session_state as state

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent


def input_example() -> None:
    """
    Set the example input for the Streamlit app.

    This function loads the example gene set and populates the `gene_set_input` and
    `gene_set_name` in the session state with the content and name of the example respectively.
    """
    logger.info("Setting the example input for the Streamlit app.")
    # Callback because that's the only way it works
    state.gene_set_input = (
        ROOT / "data" / "gene_lists" / "example_gene_list.txt"
    ).read_text()
    state.gene_set_name = "Example gene set"


def update_text_widgets() -> None:
    if "selected_file" in state and state.selected_file != "Select ...":
        file_path = ROOT / "data" / "gene_lists" / state.selected_file

        with open(file_path, "r") as file:
            file_content = file.read()

        state.gene_set_input = file_content
        state.gene_set_name = state.selected_file
