# Data Analysis README

This directory is divided into the following sub-directories, each representing a specific part of our analysis: <br>
- 'Part1_Mtb_Init_DataExploration': <br>
- 'Part2_Mtb_SVGraph': 
- 'Part3_Mtb_PG_Comparison': 
- 'Part4_BacAnno_Comparison': 
- 'Part5_Eval_SRAsm_Absences': 
- 'Part6_Ecoli_Init_DataExploration': 
- 'Part7_Ecoli_PG_Comparison': 
- 'Part8_MtbAndEcoli_NRC_Evaluation': 


The majority of analysis was performed using Python based Jupyter notebooks to process and visualize results. The directories and notebooks are ordered in the order they should be run, with earlier results in the paper coming first.

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




