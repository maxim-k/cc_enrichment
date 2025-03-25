# Transparent Gene Set Enrichment Analysis

## Prerequisites
- Have a text file with a set of newline-separated background genes in `/data/backgrounds`.
- Place gene set library files in the [Gene Matrix Transposed format](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) in `/data/libraries`.
- `alias.json` files provides human-readable names for gene set libraries and background gene sets.

## Run the Application
- To run locally, use `streamlit run code/streamlit.app`
- To run in Docker use `docker build -t cc_enrichment:latest .` and `docker run -p 8080:8080 cc_enrichment:latest`

## Overview of Application Components
### Key Features

#### Gene Set Submission
Users have the flexibility to input their gene set in two ways:
- Direct Input: Paste a newline-separated list of gene names and assign a name to the gene set.
- File Selection: Choose an existing gene list from the local `data` folder.

#### Gene Set Validation
The app performs crucial validation steps to ensure the accuracy of the analysis:
- Duplicate Check: It checks for any duplicate gene names within the input gene set.
- Background Check: It validates the input gene set against a predefined background gene set, which helps in identifying any genes that might not be present in the background set.

#### Results Display and Interaction
- Interactive Chart: The results for each gene set library are displayed in an interactive chart, allowing users to engage with the data dynamically.
- Interactive Bar Graph: An interactive bar graph provides a visual representation of the results for each gene set library.

#### Data Export Options
- Individual Library Export: Users can save the results for each gene set library in TSV and JSON formats.
- Consolidated Export: There is an option to save results for all libraries in a single TSV file.

#### Results Management
All results, along with the complementary metadata, are stored in the `results` folder and are assigned unique filenames.

### Advanced Features

#### Customization of Results Display
Users can customize the number of results to display in both the chart and the bar graph in the Streamlit app.

#### P-value Calculation Methods
The app offers various statistical methods to calculate p-values:
- Fisher's Exact Test.
- Hypergeometric Test.
- Chi-Squared Test.

#### Custom Reference Data Upload
Background Gene Set Upload: Users can upload their background gene sets as .txt files containing newline-separated gene names.
Gene Set Libraries Upload: The app allows for the upload of gene set libraries in the .gmt file format, enabling the use of custom gene sets in the analysis.
