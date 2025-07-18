import logging
from typing import Dict

import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def collect_results(results: Dict) -> str:
    """
    Concatenate all enrichment results. For each result adds a column with a library name.

    :param results: The dictionary containing enrichment results.
    """
    logger.info("Concatenating all enrichment results.")
    results_concat = []
    for library_name in results.keys():
        result = results[library_name].to_dataframe()
        result.insert(0, "Library", library_name)
        results_concat.append(result)

    return pd.concat(results_concat, ignore_index=True).to_csv(sep="\t", index=False)
