# Snakemake Pipelines README
All data processing starting from raw reads (FASTQ) was orchestrated using the Snakemake workflow system.

Each Snakemake file describes softwares & commands run (w/ specific parmeters). 
If a specific Conda environment was used for a step, a .yml file containing software/package versions will be specified.


### Overview of snakemake pipelines

#### a) Genome Assembly & Curation
1 - `Mtb.Generate.HybridAsm.PacBioRSII.smk`: Pipeline for hybrid assembly of *Mtb* isolates sequenced with both PacBio (RSII) and Illumina WGS

2 - `Mtb.Generate.HybridAsm.PacBioHiFi.smk`: Pipeline for hybrid assembly of *Mtb* isolates sequenced with both PacBio (HiFi, Sequel II) and Illumina WGS

3 - `Mtb.Generate.HybridAsm.ONT9.4.1.smk`: Pipeline for hybrid assembly of *Mtb* isolates sequenced with both Oxford Nanopore (9.4.1) and Illumina WGS

4 - `Mtb.Generate.SRAsms.smk`: Pipeline for *de novo* assembly of paired-end short reads for 151 *Mtb* isolates

5 - `Ecoli.Generate.SRAsms.smk`: Pipeline for *de novo* assembly of paired-end short reads for 50 *E. coli* isolates

#### b) *Mtb* analysis & processing

6 - `Mtb.HybridAsms.BuildPhylogeny.smk`: Pipeline for inferring a phylogeny for 151 *Mtb* isolates (Uses IQ-Tree)

7 - `Mtb.HybridAsms.Minigraph.smk`: Pipeline for generation of SV pan-genome graph of all hybrid *Mtb* assemblies

8 - `Mtb.HybridAsms.PGAnalysis.smk`: Pipeline for running genome annotation and pan-genome analysis (Hybrid *Mtb* assemblies)

9 - `Mtb.SRAsms.PGAnalysis.smk`: Pipeline for running genome annotation and pan-genome analysis (Short-read *Mtb* assemblies)

#### c) *E. coli* analysis & processing

10 - `Ecoli.HybridAsms.PGAnalysis.smk`: Pipeline for running genome annotation and pan-genome analysis (Hybrid *E. coli* assemblies)

11 - `Ecoli.SRAsms.PGAnalysis.smk`: Pipeline for running genome annotation and pan-genome analysis (Short-read *E. coli* assemblies)


### Usage
🚧 Check back soon 🚧


