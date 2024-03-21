# Manuscript Code Repository: *Analysis of the limited M. tuberculosis pangenome reveals potential pitfalls of pan-genome analysis approaches*
This repository contains all code and bioinformatics pipelines to reproduce the analysis of *Analysis of the limited M. tuberculosis pangenome reveals potential pitfalls of pan-genome analysis approaches*. <br>

All initial bioinformatics processing (such as assembly, alignment, annotation, and pan-genome estimation) was implemented using the [Snakemake](https://snakemake.github.io/) workflow engine, while downstream data processing and analysis is available in Python based Jupyter notebooks.

## Contents
- [Installation](#Installation)
- [Snakemake Pipelines](#Snakemake-pipelines)
- [Data Analysis](#Data-Analysis)
- [Results](#Results)
- [License](#License)

## Installation
All dependencies needed to reproduce this analysis can be installed via [Conda](https://docs.conda.io/en/latest/) .
🚧 Check back soon 🚧
```
# 1) Clone repository
git clone https://github.com/farhat-lab/mtb-pg-benchmarking-2024paper

# 2) Create a conda environment named 'CoreEnv_PG_V1'
cd mtb-pg-benchmarking-2024paper/

conda env create --file CondaEnvs/_____.yml -n CoreEnv_PG_V1

# 3) Activate environment (used for SnakeMake pipeline and data analysis)
conda activate CoreEnv_PG_V1
```

## Snakemake pipelines
All data processing starting from raw reads (FASTQ) was orchestrated using the Snakemake workflow system.
Refer to the Snakemake directory's [README.md](Snakemake_Pipelines/) for detailed usage of each pipeline.


## Data Analysis 

The [DataAnalysis/](https://github.com/farhat-lab/mtb-pg-benchmarking-2024paper/tree/main/Analysis) directory contains Jupyter notebooks for downstream data processing, table generation, and figure generation related to this work.

Data analysis can be run after initial assembly, annotation, and analysis pipelines have been run. (Alternatively, processed pipeline outputs can be downloaded.)

Refer to the Analysis directory's [README.md](Analysis/) for a detailed overview of this project's analysis.


## Supplemental Files
🚧 Check back soon 🚧


## Curated Datasets
All Hybrid and short-read genome assemblies (along with annotations) used in this work are organized and made easily accessible for future use.
Additionally, the (relatively) large outputs of the various pan-genome analyses are made easily available via Zenodo.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10846276.svg)](https://doi.org/10.5281/zenodo.10846276)

### Genome assemblies used
🚧 Check back soon 🚧

### Pan-genome analyses
🚧 Check back soon 🚧

## License
This repository is distributed under the [MIT license terms](LICENSE).
