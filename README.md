# Manuscript Code Repository: *Analysis of the limited M. tuberculosis pangenome reveals potential pitfalls of pan-genome analysis approaches*
This repository contains all code and bioinformatics pipelines to reproduce the analysis of *Analysis of the limited M. tuberculosis pangenome reveals potential pitfalls of pan-genome analysis approaches*. <br>

All initial bioinformatics processing (such as assembly, alignment, annotation, and pan-genome estimation) is implemented using the  [Snakemake](https://snakemake.github.io/) workflow system, while downstream data processing and analysis is available in Jupyter notebooks.

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
```

## Snakemake pipelines
All data processing starting from raw reads (FASTQ) was orchestrated using the Snakemake workflow system.
Refer to the Snakemake directory's [README.md](Snakemake_Pipelines/) for detailed usage of each pipeline.


## Data Analysis 

The [DataAnalysis/](https://github.com/farhat-lab/mtb-illumina-wgs-evaluation/tree/main/DataAnalysis) directory contains Jupyter notebooks for downstream data processing, table generation, and figure generation related to this work.

Data analysis can be run after initial assembly, annotation, and analysis pipelines have been run.

### Overview of phases of data analysis:
1) Summarize isolate and assembly (Hybrid & SR) info for Mtb dataset
2) Comparing assembly stats for Hybrid & SR assemblies
3) Phylogeny visualization & pairwise genome similarity comparison
4) Processing & analysis of SV Mtb pan-genome graph
5) Analysis of *Mtb* pan-genome estimates
6) Comparison of annotation pipelines (`Bakta`, `PGAP`)
7) Benchmarking pan-genome estimates based on SR assemblies (Leveraging available hybrid assemblies as ground truth)
8) Analysis of NRC applied to *Mtb* (Adjusting for DNA seq redundancy in pan-genome estimates)
9) Organization of *E. coli* dataset (From [Shaw et. al. 2021](https://www.science.org/doi/10.1126/sciadv.abe3868))
10) Analysis of *E. coli* pan-genome estimates
11) Analysis of NRC applied to *E. coli* (Adjusting for DNA seq redundancy in pan-genome estimates)


## Supplemental Files
🚧 Check back soon 🚧


## Curated Datasets
All Hybrid and short-read genome assemblies (w/ annotations) used in this work are organized and made easily accessible for future use.
Additionally, the (relatively) large outputs of the various pan-genome analyses are made easily available via Zenodo & Harvard Dataverse.

### Genome assemblies used
🚧 Check back soon 🚧

### Pan-genome analyses
🚧 Check back soon 🚧

## License
This repository is distributed under the [MIT license terms](LICENSE).
