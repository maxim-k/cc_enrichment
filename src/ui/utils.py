import base64
import json
import logging
import re
from pathlib import Path
from pprint import pformat
from typing import Dict

import streamlit as st

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parent.parent.parent


def update_aliases(directory: str, alias_file: str = "alias.json") -> Dict[str, str]:
    """
    Update the aliases file with the current set of files in the specified directory.

    This function scans a specified directory for files, keeping an 'alias.json' file that contains
    a mapping of simplified names to actual file names. This allows for easier reference to the files
    in that directory. If new files are added to the directory, they are added to the 'alias.json'. If files
    are missing, they are removed from 'alias.json'.

    :param directory: The directory to scan for files.
    :param alias_file: The name of the alias file, defaults to 'alias.json'
    :return: A dictionary containing the aliases, with keys being the simplified names and values being the actual file names.
    """
    logger.info(f"Updating aliases for directory: {directory}")
    aliases_path = ROOT / "data" / directory / alias_file

    if Path(aliases_path).is_file():
        try:
            with open(aliases_path, "r") as file:
                alias = json.load(file)
        except (FileNotFoundError, json.JSONDecodeError):
            logger.warning(f"Failed to load aliases from {aliases_path}")
            st.warning(f"Failed to load aliases from {aliases_path}")

    files = [f for f in (ROOT / "data" / directory).iterdir() if f.is_file()]

    # Remove 'alias.json' from the list of files
    if Path(aliases_path) in files:
        files.remove(Path(aliases_path))

    # Add a record if a file is not in aliases
    for file in files:
        if file.name not in alias.values():
            alias[file.stem] = file.name

    # Delete a record from aliases if there's no corresponding file
    aliases_keys_to_delete = [
        key for key in alias if alias[key] not in [file.name for file in files]
    ]

    for key in aliases_keys_to_delete:
        del alias[key]

    with open(aliases_path, "w") as file:
        json.dump(alias, file, indent=4)
    logger.info(f"{directory}/alias.json\n{pformat(alias)}")
    return alias


def download_link(val: str, filename: str, extension: str) -> str:
    """
    Create a download link for a file with the given content, filename, and extension.

    This function generates a download link that, when clicked, will download the file with the given
    content, filename, and extension. The content is encoded in base64 and the link is created using an
    HTML 'a' tag with a 'download' attribute.

    :param val: The content of the file to be downloaded.
    :param filename: The name of the file, without the extension.
    :param extension: The file extension (e.g., 'tsv', 'json').
    :return: An HTML string containing the download link.
    """
    logger.info(f"Creating download link for file: {filename}.{extension}")
    b64 = base64.b64encode(val.encode("utf-8"))
    return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{filename}.{extension}">{extension}</a>'


def sanitize_id(raw: str) -> str:
    """
    Convert raw label into a valid DOT node ID: replace non-alphanumeric with underscores.
    Collapse multiple underscores and strip leading/trailing underscores.
    """
    # replace non-word characters with underscore
    s = re.sub(r"\W+", "_", raw)
    # collapse underscores
    s = re.sub(r"_+", "_", s)
    return s.strip("_")
