# **🚧 Repo under construction 🚧**

# Genomic and Epidemiologic Factors Unify Distinct Esophageal Adenocarcinoma 

This repository provides the code, data, and supplementary materials for the study:

**"Genomic and Epidemiologic Factors Unify Distinct Esophageal Adenocarcinoma" [**Link to the paper or DOI] Published in *[Add Journal Name]*, [Year].

## Authors and Affiliations {#authors-and-affiliations}

**Authors:**\
- Shahriar A. Zamani¹²³*, Lianlian Wu¹*, Emily L. Black¹, Alex Bartram¹, Alvin W. T Ng¹⁴, Maria Secrier⁵, Daniel Jacobson¹, Ginny Devonshire¹, Nicola Grehan¹, Barbara Nutzinger¹, Adam Freeman¹, Ahmad Miremadi⁶, Maria O’Donovan⁶, Alex M. Frankell¹, Sarah Killcoyne¹, Oesophageal Cancer Clinical and Molecular Stratification (OCCAMS) Consortium, Helen G. Coleman⁷, Rebecca C. Fitzgerald¹#

**Affiliations:**\
1. Early Cancer Institute, Department of Oncology, University of Cambridge, Cambridge, England, United Kingdom.\
2. Cancer Prevention Fellowship Program, Division of Cancer Prevention, National Cancer Institute, National Institutes of Health, Rockville, Maryland, United States.\
3. Division of Cancer Epidemiology and Genetics, National Cancer Institute, Bethesda, Maryland, United States.\
4. Cancer Research UK Cambridge Institute, University of Cambridge, Cambridge, England, United Kingdom.\
5. UCL Genetics Institute, Department of Genetics, Evolution and Environment, University College London, London, England, United Kingdom.\
6. Department of Pathology, Cambridge University Hospital NHS Trust, Cambridge, England, United Kingdom.\
7. Cancer Epidemiology Research Group, Centre for Public Health, Queen's University Belfast, Belfast, Northern Ireland, United Kingdom.

**Notes:**\
- \* Dual first-authors\
- \# Corresponding author

## Overview

This study investigates the genomic and epidemiologic characteristics of esophageal adenocarcinoma (EAC) with a focus on whether Barrett's Esophagus (BE) is a necessary precursor. By analyzing data from a large, prospective cohort, the findings suggest that BE is the precursor lesion to EAC, though it may be obscured in advanced cases.

## Table of Contents

-   [Installation](#installation)
-   [Usage](#usage)
-   [Repository Structure](#repository-structure)
-   [Citation](#citation)
-   [License](#license)
-   [Acknowledgments](#acknowledgments)

## Installation {#installation}

Clone this repository and install the necessary Python dependencies.

``` bash
git clone https://github.com/yourusername/repo-name.git
cd repo-name
pip install -r requirements.txt
```

[Add additional setup steps if necessary.]

## Usage {#usage}

This repository provides tools for: 1. **Data Analysis:** Reproducing the logistic regression models and multiregional mutational lineage tracing. 2. **Data Visualization:** Generating plots and phylogenetic trees included in the manuscript. 3. **Supplementary Material:** Accessing and exploring supplementary figures and tables.

### Example Commands

To run the logistic regression analysis:

``` bash
python src/logistic_regression.py --input data/cohort_data.csv --output results/analysis_results.csv
```

To generate a phylogenetic tree:

``` bash
python src/phylogenetics.py --input data/sequencing_data.csv --output figures/phylogenetic_tree.png
```

## Repository Structure {#repository-structure}

```         
repo-name/
├── data/               # Example datasets (de-identified)
├── src/                # Scripts for analysis and visualization
├── figures/            # Generated figures
├── results/            # Output from analyses
├── notebooks/          # Jupyter notebooks for detailed exploratory analyses
├── requirements.txt    # Python dependencies
└── README.md           # Project description
```

## Citation {#citation}

If you use this repository, please cite our paper:

```         
@article{zamani2024,
  author    = {Shahriar A. Zamani and Lianlian Wu and others},
  title     = {Genomic and Epidemiologic Factors Unify Distinct Esophageal Adenocarcinoma},
  journal   = {Add Journal Name},
  year      = {2024},
  volume    = {XX},
  pages     = {XX--XX},
  doi       = {10.xxxx/xxxx},
}
```

## License {#license}

This project is licensed under the [MIT License](LICENSE).

## Acknowledgments {#acknowledgments}

This study was supported by Cancer Research UK, the Medical Research Council, and other funding bodies listed in the manuscript. Special thanks to the OCCAMS Consortium, study participants, and collaborators who made this work possible.
