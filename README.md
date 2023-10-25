# CodeOcean Gene Set Enrichment Analysis

## Prerequisites
- Have a text file with a set of newline-separated background genes in `/data/backgrounds`.
- Place gene set library files in the [Gene Matrix Transposed format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) in `/data/libraries`.
- `alias.json` files provides human-readable names for gene set libraries and background gene sets.

## Run the App
- To run locally, use `streamlit run code/stramlit.app`
- To run on CodeOcean, click the Streamlit icon on the Cloud Workstation dashboard.