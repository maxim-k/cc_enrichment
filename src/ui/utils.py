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
    Scan the given directory and update alias.json, supporting both the old dict format
    and the new list-of-dicts format. Returns a mapping of alias -> filename for active entries.
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Updating aliases for directory: {directory}")
    aliases_path = ROOT / "data" / directory / alias_file

    alias_entries = []
    is_list_format = False

    # Load existing aliases
    if aliases_path.is_file():
        try:
            with open(aliases_path, "r") as f:
                data = json.load(f)
            if isinstance(data, list):
                # New format: list of {name, file, active}
                is_list_format = True
                alias_entries = data
            elif isinstance(data, dict):
                # Old format: simple mapping
                alias_entries = [
                    {"name": name, "file": filename, "active": True}
                    for name, filename in data.items()
                ]
            else:
                logger.warning(f"Unrecognized alias.json format: {type(data)}")
        except (FileNotFoundError, json.JSONDecodeError):
            logger.warning(f"Failed to load aliases from {aliases_path}")
            st.warning(f"Failed to load aliases from {aliases_path}")

    # Collect actual files in the directory (excluding the alias file)
    dir_path = ROOT / "data" / directory
    files = [
        f for f in dir_path.iterdir()
        if f.is_file() and f.name != alias_file
    ]

    # Track existing filenames in entries
    existing_files = {entry["file"] for entry in alias_entries}

    # Add new files as active
    for file in files:
        if file.name not in existing_files:
            alias_entries.append({
                "name": file.stem,
                "file": file.name,
                "active": True
            })

    # Remove entries whose files no longer exist
    current_files = {f.name for f in files}
    alias_entries = [
        entry for entry in alias_entries
        if entry["file"] in current_files
    ]

    # Build result mapping for active entries only
    alias_mapping: Dict[str, str] = {
        entry["name"]: entry["file"]
        for entry in alias_entries
        if entry.get("active", False)
    }

    # Write back to alias.json in the original format
    with open(aliases_path, "w") as f:
        if is_list_format:
            json.dump(alias_entries, f, indent=4)
        else:
            json.dump(alias_mapping, f, indent=4)

    logger.info(f"{directory}/alias.json updated: {alias_mapping}")
    return alias_mapping


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
