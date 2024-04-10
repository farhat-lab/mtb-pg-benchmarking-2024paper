# 8.Mtb.HybridAsms.PGAnalysis.smk
### Snakemake pipeline for performing pan-genome analysis of all Mtb hybrid assemblies (N = 151)


### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]

H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]
H37rv_GBK_PATH = config["H37rv_GBK_PATH"]


# Define PATH of main OUTPUT directory
output_Dir = config["output_dir"]




# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')


# Create a python list of Sample IDs
input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )

SampleID_to_Asm_FA_Dict = dict(input_DataInfo_DF[['SampleID', 'Genome_ASM_PATH']].values)

print("List of input sampleIDs:", len(input_All_SampleIDs), input_All_SampleIDs)



# Define directory with PGAP Annotated Assemblies
PGAP_GenomeAnnoDir = "/n/data1/hms/dbmi/farhat/mm774/Projects/230121_PGAP_AnnoByTBPortals_V1/genomes"

PGAP_SRAsm_GenomeAnnoDir = "/n/data1/hms/dbmi/farhat/mm774/Projects/230313_PGAP_AnnoByTBPortals_SRAsms_V1/genomes"





rule all:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.fna", sampleID = input_All_SampleIDs),
        output_Dir + "/AsmAnalysis/H37Rv/GenomeAnnotation/H37Rv_Asm_Bakta/H37Rv.Bakta.gff3",

        expand(output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_LRAsm/{sampleID}/short_summary.specific.corynebacteriales_odb10.{sampleID}.txt", sampleID = input_All_SampleIDs), #["MFS-56"]),
        output_Dir + "/AsmAnalysis/H37Rv/Assembly_QC/Busco_Mtb_LRAsm/H37Rv/short_summary.specific.corynebacteriales_odb10.H37Rv.txt",
        
        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/pangenome.ContentSummary.txt",

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/Panstripe_PhyloAnalysis/Panstripe.PanstripeSummary.tsv"


        # PG Analysis based on PGAP Annotations
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = ["M0011368_9"]),
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.gff", sampleID = input_All_SampleIDs), #input_All_SampleIDs

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.ContentSummary.txt",


        # Run panqc (Nucleotide Redundancy & CDS Annotation Discrepancy Adjustment/Analysis)
        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",


        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",




##################################################
############ Genome Assembly QC ##################
##################################################

# Run BUSCO for general QC and completeness evaluation

rule busco_LRAsm:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        #Rv_ShortSum_TXT = output_Dir + "/AsmAnalysis/H37Rv/Assembly_QC/Busco_Mtb_LRAsm/H37Rv/short_summary.specific.corynebacteriales_odb10.H37Rv.txt",
    output:
        ShortSum_TXT = output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_LRAsm/{sampleID}/short_summary.specific.corynebacteriales_odb10.{sampleID}.txt",
    conda:
        "CondaEnvs/Busco_v5_4_4.nobuilds.yml"
    params:
        mode = "genome",
        name = '{sampleID}',
        out_dir = output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_LRAsm/",
        tmp_dir = output_Dir + "/Busco_Download_Tmp/",
        lineage="corynebacteriales_odb10",
    shell:
        #"mkdir {params.tmp_dir} \n"
        "mkdir -p {params.out_dir} \n"
        "cd {params.out_dir} \n"
        "busco -i {input.i_Assembly_FA} -o {wildcards.sampleID} --out_path {params.out_dir} --mode {params.mode} "
        " --lineage_dataset {params.lineage} --download_path {params.tmp_dir} --force --offline "

rule busco_H37Rv:
    input:
        H37Rv_Assembly_FA = refGenome_FA_PATH
    output:
        Rv_ShortSum_TXT = output_Dir + "/AsmAnalysis/H37Rv/Assembly_QC/Busco_Mtb_LRAsm/H37Rv/short_summary.specific.corynebacteriales_odb10.H37Rv.txt",
    conda:
        "CondaEnvs/Busco_v5_4_4.nobuilds.yml"
    params:
        mode = "genome",
        name = 'H37Rv',
        out_dir = output_Dir + "/AsmAnalysis/H37Rv/Assembly_QC/Busco_Mtb_LRAsm/",
        tmp_dir = output_Dir + "/Busco_Download_Tmp/",
        lineage="corynebacteriales_odb10",
    shell:
        "mkdir -p {params.tmp_dir} \n"
        "busco -i {input} -o {params.name} --out_path {params.out_dir} --mode {params.mode} "
        " --lineage_dataset {params.lineage} --download_path {params.tmp_dir} --force "


#################################################
############ Genome Annotation ##################
#################################################


rule Bakta_Anno_WiH37Rv:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        i_H37Rv_GBK = H37rv_GBK_PATH,
    output:
        Asm_Bakta_Anno_FAA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.faa",
        Asm_Bakta_Anno_FNA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.fna",
        Asm_Bakta_Anno_GBFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gbff",
        Asm_Bakta_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3",
    conda:
        "CondaEnvs/Bakta_1_6_1.nobuilds.yml" 
    threads: 8
    params:
        Bakta_OutputDir_PATH = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/",
        Bakta_DB_Dir = "/n/data1/hms/dbmi/farhat/mm774/References/Bakta_DB_v4/db",
        Bakta_OutPrefix = "{sampleID}.Bakta"
    shell:
        "bakta --db {params.Bakta_DB_Dir} --verbose --output {params.Bakta_OutputDir_PATH} "
        " --prefix {params.Bakta_OutPrefix} --locus-tag {wildcards.sampleID}  "
        " --proteins {input.i_H37Rv_GBK}"
        " --threads {threads} {input.i_Assembly_FA} "




###############################################


# https://github.com/gtonkinhill/panaroo/issues/68
rule CleanUp_Bakta_GFF:
    input:
        Asm_Bakta_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3",
        Asm_Bakta_Anno_FNA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.fna",
    output:
        Asm_Bakta_Anno_IDsFixed_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.IDs.Fixed.gff3",
        Asm_Bakta_Anno_IDsFixed_FNA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.IDs.Fixed.fna",
    shell:
        'grep -v "ID[123]=" {input.Asm_Bakta_Anno_GFF} > {output.Asm_Bakta_Anno_IDsFixed_GFF} \n'
        "cp  {input.Asm_Bakta_Anno_FNA} {output.Asm_Bakta_Anno_IDsFixed_FNA} "


rule Anno_H37Rv_WiBakta:
    input:
        H37rv_Ref_fa = refGenome_FA_PATH,
        i_H37Rv_GBK = H37rv_GBK_PATH,
    output:
        Asm_Bakta_Anno_FAA = output_Dir + "/AsmAnalysis/H37Rv/GenomeAnnotation/H37Rv_Asm_Bakta/H37Rv.Bakta.faa",
        Asm_Bakta_Anno_FNA = output_Dir + "/AsmAnalysis/H37Rv/GenomeAnnotation/H37Rv_Asm_Bakta/H37Rv.Bakta.fna",
        Asm_Bakta_Anno_GBFF = output_Dir + "/AsmAnalysis/H37Rv/GenomeAnnotation/H37Rv_Asm_Bakta/H37Rv.Bakta.gbff",
        Asm_Bakta_Anno_GFF = output_Dir + "/AsmAnalysis/H37Rv/GenomeAnnotation/H37Rv_Asm_Bakta/H37Rv.Bakta.gff3",
    conda:
        "CondaEnvs/Bakta_1_6_1.nobuilds.yml"
    threads: 8
    params:
        Bakta_OutputDir_PATH = output_Dir + "/AsmAnalysis/H37Rv/GenomeAnnotation/H37Rv_Asm_Bakta/",
        Bakta_DB_Dir = "/n/data1/hms/dbmi/farhat/mm774/References/Bakta_DB_v4/db",
        Bakta_OutPrefix = "H37Rv.Bakta"
    shell:
        "bakta --db {params.Bakta_DB_Dir} --verbose --output {params.Bakta_OutputDir_PATH} "
        " --prefix {params.Bakta_OutPrefix} --locus-tag H37Rv  "
        " --proteins {input.i_H37Rv_GBK}"
        " --threads {threads} {input.H37rv_Ref_fa} "







##################################################
############ Pangenome Analysis ##################
##################################################


rule Panaroo_StdParams_Strict_AllIsolates_WiH37Rv:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),
       H37rv_Ref_GFF = "/n/data1/hms/dbmi/farhat/mm774/References/GCF_000195955.2_ASM19595v2_genomic.ForPanaroo.gff3",            
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiH37Rv/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiH37Rv/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiH37Rv/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiH37Rv/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes"



rule Panaroo_StdParams_Strict_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes"


rule Panaroo_StdParams_Moderate_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes"


rule Panaroo_StdParams_Sensitive_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes"




##### Run Panaroo but with the '--merge_paralogs' option. 

rule Panaroo_StdParams_Strict_MergeParalogs_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Moderate_MergeParalogs_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Sensitive_MergeParalogs_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",

    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes --merge_paralogs"





rule Roary_DefaultParams_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        'TMPDIR={params.Roary_TmpDir}  \n'
        'mkdir -p {params.Roary_TmpDir}  \n'
        'echo $TMPDIR \n'
        #"roary -r -p {threads} -f {params.Roary_OutputDir} {input} "
        "roary -e -n -r -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e --mafft


rule Roary_NoSplitParalogs_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 1
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        #'TMPDIR={params.Roary_TmpDir}  \n'
        #'mkdir -p {params.Roary_TmpDir}  \n'
        #'echo $TMPDIR \n'
        #"roary -r -s -p {threads} -f {params.Roary_OutputDir} {input} " # -s means 'dont split paralogs'
        "roary -e -n -v -r -s -p {threads} -f {params.Roary_OutputDir} {input} " # -s means 'dont split paralogs'
        # -e --mafft



rule Roary_NoSplitParalogs_I90_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        'TMPDIR={params.Roary_TmpDir}  \n'
        'mkdir -p {params.Roary_TmpDir}  \n'
        'echo $TMPDIR \n'
        #"roary -r -i 90 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        "roary -e -n -r -i 90 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -s means 'dont split paralogs'
        # -i 90 sets the minimumblastp identity to 90%
        # -e --mafft


rule Roary_NoSplitParalogs_I80_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        #'TMPDIR={params.Roary_TmpDir}  \n'
        #'mkdir -p {params.Roary_TmpDir}  \n'
        #'echo $TMPDIR \n'
        #"roary -r -i 80 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        "roary -e -n -r -i 80 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -s means 'dont split paralogs'
        # -i 80 sets the minimumblastp identity to 80%
        # -e --mafft
 


rule create_InputPATHInfo_Bakta_GBFF:
    input:
        Asm_Bakta_Anno_GBFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gbff",
    output:
        Bakta_GBFF_InputPATH_Info_TSV = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gbff.InputInfo.tsv",
    shell:
        'echo "{wildcards.sampleID}\t{input.Asm_Bakta_Anno_GBFF}" > {output.Bakta_GBFF_InputPATH_Info_TSV}'


rule mergeInputGBFF_PATHs_To_TSV_ForPpanggolin_AllIsolates:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gbff.InputInfo.tsv", sampleID = input_All_SampleIDs),       
    output:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/Ppanggolin_Preprocessing/AllGenomes.GBFF.InputPATHs.tsv",
    shell:
        "cat {input} > {output.All_Input_GBFF_PATHs_TSV}"


rule run_Ppanggolin_DefaultParam_AllIsolates:
    input:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/Ppanggolin_Preprocessing/AllGenomes.GBFF.InputPATHs.tsv",
    output:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/pangenome.h5",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 8
    params:
        Ppanggolin_OutputDir = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/",
    shell:
        "ppanggolin workflow -f --cpu {threads} --output {params.Ppanggolin_OutputDir} --anno {input.All_Input_GBFF_PATHs_TSV} "

rule get_PangenomeContent_Ppanggolin_AllIsolates:
    input:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/pangenome.h5",
    output:
        Pangenome_ContentSummary_TXT = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/pangenome.ContentSummary.txt",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 1
    shell:
        "ppanggolin info -p {input.Pangenome_H5} --content > {output.Pangenome_ContentSummary_TXT}"


rule get_PG_GeneSeqs_Ppanggolin_AllIsolates:
    input:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/pangenome.h5",
    output:
        Pangenome_GeneSeqs_FA = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/all_nucleotide_families.fasta",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 1
    params:
        Ppanggolin_OutputDir = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates/",
    shell:
        "ppanggolin fasta -p {input.Pangenome_H5} --output {params.Ppanggolin_OutputDir} --gene_families all --force"


##################################################






##################################################################################
##################### Run PGAP Based Pangenome Analysis  #########################
##################################################################################

# /n/data1/hms/dbmi/farhat/mm774/Projects/230121_PGAP_AnnoByTBPortals_V1/genomes/R37765/output


rule CP_PGAP_AnnoFiles:
    input:
        i_Anno_FNA = PGAP_GenomeAnnoDir + "/{sampleID}/output/annot.fna",      
        i_Anno_GFF = PGAP_GenomeAnnoDir + "/{sampleID}/output/annot.gff", 
        i_Anno_WiGenomicFASTA_GFF = PGAP_GenomeAnnoDir + "/{sampleID}/output/annot_with_genomic_fasta.gff",        
        i_Anno_GBK = PGAP_GenomeAnnoDir + "/{sampleID}/output/annot.gbk",  
        i_Anno_Genome_FA = PGAP_GenomeAnnoDir + "/{sampleID}/output/annot.fna",  
    output:
        Renamed_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.gff", 
        Renamed_Genome_FA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.Genome.fasta",  

    shell: # Rename the GFF & FASTA file to have the SampleID as the contig name
        "awk -F'\\t' '/^#/ {{print; next}} {{OFS=\"\\t\"; $1=\"{wildcards.sampleID}\"; print}}' {input.i_Anno_GFF} >  {output.Renamed_Anno_GFF} \n"
        " bioawk -c fastx '{{ print \">{wildcards.sampleID}\" \"\\n\" $seq }}' {input.i_Anno_Genome_FA} > {output.Renamed_Genome_FA}"

rule reformat_PGAP_GFF:
    input:
        Renamed_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.gff", 
        Renamed_Genome_FA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.Genome.fasta",  
    output:
        Updated_PGAP_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", 
    shell:
        "convert_refseq_to_prokka_gff.py -g {input.Renamed_Anno_GFF} -f {input.Renamed_Genome_FA} -o {output.Updated_PGAP_GFF} "


# Run Pangenome analysis tools w/ PGAP Anno Asms

#  PGAP - Panaroo
rule Panaroo_StdParams_Strict_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
       #expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 4
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes"



rule Panaroo_StdParams_Moderate_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes"



rule Panaroo_StdParams_Sensitive_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes"



# PGAP - Panaroo - Wi Merge Paralogs

##### Run Panaroo but with the '--merge_paralogs' option. 

rule Panaroo_StdParams_Strict_MergeParalogs_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes --merge_paralogs"



rule Roary_DefaultParams_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates_WiPGAPAnno_V1/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary  -v -r -p {threads} -f {params.Roary_OutputDir} {input} " 
        # -e -n These flags create a `pan_genome_reference.fa` file


rule Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary -v -r -s -p {threads} -f {params.Roary_OutputDir} {input} " # -s means 'dont split paralogs'
        # -e -n These flags create a `pan_genome_reference.fa` file


rule Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "export R_DEFAULT_SAVE_WORKSPACE=FALSE \n"
        'rm -r {params.Roary_OutputDir} \n'
        "roary -v -r -i 90 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e -n These flags create a `pan_genome_reference.fa` file

        # -s means 'dont split paralogs'
        # -i 90 sets the minimumblastp identity to 90%
        # -e --mafft


rule Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "export R_DEFAULT_SAVE_WORKSPACE=FALSE \n"
        'rm -r {params.Roary_OutputDir} \n'
        "roary -v -r -i 80 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e -n These flags create a `pan_genome_reference.fa` file

        # -s means 'dont split paralogs'
        # -i 80 sets the minimumblastp identity to 80%
        # -e --mafft
 


rule create_InputPATHInfo_PGAPAnno_GFF:
    input:
        Asm_PGAP_Anno_GBFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff",
    output:
        PGAP_GFF_InputPATH_Info_TSV = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.gff.InputInfo.tsv",
    shell:
        'echo "{wildcards.sampleID}\t{input.Asm_PGAP_Anno_GBFF}" > {output.PGAP_GFF_InputPATH_Info_TSV}'


rule mergeInputGBFF_PATHs_To_TSV_ForPpanggolin_AllIsolates_WiPGAPAnno:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_PGAP_V1/{sampleID}.PGAP.gff.InputInfo.tsv", sampleID = input_All_SampleIDs),       
    output:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/Ppanggolin_Preprocessing_PGAP/AllGenomes.GFF.InputPATHs.tsv",
    shell:
        "cat {input} > {output.All_Input_GBFF_PATHs_TSV}"


rule run_Ppanggolin_DefaultParam_AllIsolates_WiPGAPAnno:
    input:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/Ppanggolin_Preprocessing_PGAP/AllGenomes.GFF.InputPATHs.tsv",
    output:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.h5",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 8
    params:
        Ppanggolin_OutputDir = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "ppanggolin workflow -f --cpu {threads} --output {params.Ppanggolin_OutputDir} --anno {input.All_Input_GBFF_PATHs_TSV} "

rule get_PangenomeContent_Ppanggolin_AllIsolates_WiPGAPAnno:
    input:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.h5",
    output:
        Pangenome_ContentSummary_TXT = output_Dir + "/PanGenome_Analysis/Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.ContentSummary.txt",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 1
    shell:
        "ppanggolin info -p {input.Pangenome_H5} --content > {output.Pangenome_ContentSummary_TXT}"




##### Run Panstripe Phylo Analysis w/ GAIN + LOSS measure #####

# Panstripe rule in progress here!!!!

rule RunPanstripe_Panaroo_Strict_MergeParalogs:
    input:
       PA_Rtab = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.Rtab",
       in_tree = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/IQtree.10AmbThresh.PLCMask.MidRoot.WiNodeNames.newick",
    output:
       PA_NoBakta_Rtab = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.Renamed.Rtab",
       Panstripe_Summ_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/Panstripe_PhyloAnalysis/Panstripe.PanstripeSummary.tsv"
    threads: 1
    params:
        Panstripe_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/Panstripe_PhyloAnalysis",
    shell:
        'rm -r {params.Panstripe_OutputDir} \n'
        "~/conda3/bin/conda deactivate \n"
        "module load gcc/6.2.0 R/4.1.1 \n"
        'export PATH="/n/data1/hms/dbmi/farhat/mm774/Snakemake_Pipelines/Mtb-WGA-SMK/Scripts/Panstripe_Scripts:$PATH" \n'
        ""
        " sed 's/\.Bakta//g' {input.PA_Rtab} > {output.PA_NoBakta_Rtab} \n"
        "RunPhyloAnalysis.WiPanstripe.V1.R -p {input.in_tree} -r {output.PA_NoBakta_Rtab} -o {params.Panstripe_OutputDir} -a Panstripe"



#######################################################












rule Roary_DefaultParams_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "panqc ava -k {params.k} -i {input} -o {output} "

rule Roary_NoSplitParalogs_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "panqc ava -k {params.k} -i {input} -o {output} "


rule Roary_NoSplitParalogs_I90_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "panqc ava -k {params.k} -i {input} -o {output} "



rule Roary_NoSplitParalogs_I80_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "panqc ava -k {params.k} -i {input} -o {output} "





###### Run GeneSeqCheck step for all Panaroo & Roary analyses ######

rule create_SampleToAsmFA_TSV:
    output:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
    run:
        # Output as tsv
        input_DataInfo_DF[['SampleID', 'Genome_ASM_PATH']].to_csv(output.asmfa_tsv, sep="\t", index = False)

rule create_SampleToAsmFA_LRandSR_TSV:
    output:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.LRandSR.FA.PATHs.tsv",
    run:
        input_DataInfo_LRandSR_DF = input_DataInfo_DF.copy()
        input_DataInfo_LRandSR_DF.columns = ['SampleID', 'Dataset_Tag',
                                             'Genome_LR_ASM_PATH', 'Genome_SR_ASM_PATH']
        # Output as tsv
        input_DataInfo_DF[['SampleID', 'Genome_LR_ASM_PATH', 'Genome_SR_ASM_PATH']].to_csv(output.asmfa_tsv, sep="\t", index = False)





rule Panaroo_Strict_MergeParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Panaroo_Moderate_MergeParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Panaroo_Sensitive_MergeParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Roary_DefaultParams_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",  
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "

rule Roary_NoSplitParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Roary_NoSplitParalogs_I90_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Roary_NoSplitParalogs_I80_panqc_AsmGeneSeqCheck: 
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



#####################################################################


##### 

rule Panaroo_Strict_MergeParalogs_panqc_NucSimCluster:
    input:
        PG_AvA_MaxJC_TSV          = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    output:
        PresAbs_ASC_NSC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AdjBy.ASC.NSC.Rtab",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/NSC_08.Clusters.tsv",
    params:
        Ksim_ClusterThresh = 0.8,
        runtime = '0-00:6:00',
        partition = 'short'
    resources:
        mem_mb = 800,
        runtime_s = 360
    shell:
        "time panqc nscluster " 
        " --in_ava_tsv {input.PG_AvA_MaxJC_TSV} "
        " --in_gene_matrix_tsv {input.PresAbs_AsmSeqCheck_TSV} "
        " --min_ksim {params.Ksim_ClusterThresh} "
        " --out_nsc_gene_matrix {output.PresAbs_ASC_NSC_TSV} "
        " -c {output.Clusters_TSV} "


rule Panaroo_Sensitive_MergeParalogs_panqc_NucSimCluster:
    input:
        PG_AvA_MaxJC_TSV          = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    output:
        PresAbs_ASC_NSC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AdjBy.ASC.NSC.Rtab",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/NSC_08.Clusters.tsv",
    params:
        Ksim_ClusterThresh = 0.8,
        runtime = '0-00:6:00',
        partition = 'short'
    resources:
        mem_mb = 800,
        runtime_s = 360
    shell:
        "time panqc nscluster " 
        " --in_ava_tsv {input.PG_AvA_MaxJC_TSV} "
        " --in_gene_matrix_tsv {input.PresAbs_AsmSeqCheck_TSV} "
        " --min_ksim {params.Ksim_ClusterThresh} "
        " --out_nsc_gene_matrix {output.PresAbs_ASC_NSC_TSV} "
        " -c {output.Clusters_TSV} "


##### Run end-to-end panQC-NRC pipeline #####





rule Panaroo_Strict_MergeParalogs_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/panqc_NRC_OutDir",
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "


rule Panaroo_Moderate_MergeParalogs_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/panqc_NRC_OutDir",
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "



rule Panaroo_Sensitive_MergeParalogs_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/panqc_NRC_OutDir",
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "



rule Roary_DefaultParams_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir", 
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "



rule Roary_NoSplitParalogs_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/panqc_NRC_OutDir", 
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "



rule Roary_NoSplitParalogs_I90_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/panqc_NRC_OutDir", 
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "


rule Roary_NoSplitParalogs_I80_panqc_NRC:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
        PresAbs_NRC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/panqc_NRC_OutDir/gene_presence_absence.NRCUpdated.tsv",
        PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/gene_presence_absence.AsmSeqCheck.tsv",
        ASM_Stats_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/panqc_NRC_OutDir/Step1_AsmSeqCheck/AsmSeqCheck.Stats.tsv",
        Clusters_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/NSC.ClusterInfo.tsv",
        PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/panqc_NRC_OutDir/Step2_SeqClustering/AllVsAll.KmerSimilarity.tsv",
    threads: 1
    params:
        NRC_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/panqc_NRC_OutDir", 
        Ksim_ClusterThresh = 0.8,
        k = 31, # k-mer size in bp
        runtime = '0-00:10:00',
        partition = 'short',
    resources:
        mem_mb = 1000,
        runtime_s = 600
    shell:
        "time panqc nrc "
        " -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {params.NRC_OutputDir} "








#########################################################
############ LR Asm FastA ANI Analysis ##################
#########################################################

rule FastANI_I3_PP_Assembly: 
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        H37rv_Ref_fa = refGenome_FA_PATH,     
    output:
        output_Dir + "/AsmAnalysis/{sampleID}/FastANI/FastANI_AsmToH37Rv/{sampleID}.AsmToH37Rv.FastANI.txt"
    conda:
        "CondaEnvs/fastani_1_3_0_Conda.yml"
    shell:
        "fastANI -q {input.i_Assembly_FA} -r {input.H37rv_Ref_fa} -o {output}"


rule Save_LR_Asm_PATH_To_TXT: 
    input:
        Asm_Bakta_Anno_FNA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.fna",
    output:
        output_Dir + "/AsmAnalysis/{sampleID}/FastANI/{sampleID}.PathToFASTA.LRAsm.txt"
    shell:
        'echo "{input}" > {output}' 


rule Save_H37Rv_FA_PATH_To_TXT: 
    input:
        H37rv_Ref_fa = refGenome_FA_PATH,     
    output:
        output_Dir + "/FastANI/H37Rv.PathToFASTA.txt"
    shell:
        'echo "{input}" > {output}' 


rule combine_All_LR_Asm_FA_PATHs:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/FastANI/{sampleID}.PathToFASTA.LRAsm.txt", sampleID = input_All_SampleIDs),           
    output:
        output_Dir + "/FastANI/FastANI_LRAsms/LRAsms.PathToFASTAs.txt"
    shell:
        "cat {input} > {output}"


rule FastANI_LRAsm_AllVsAll: 
    input:
        output_Dir + "/FastANI/FastANI_LRAsms/LRAsms.PathToFASTAs.txt"
    output:
        output_Dir + "/FastANI/FastANI_LRAsm/FastANI.AllVsAll.LRAsm.txt"
    conda:
        "CondaEnvs/fastani_1_3_0_Conda.yml"
    threads: 8
    shell:
        "fastANI -t {threads} --ql {input} --rl {input} -o {output}"

#########################################################





#########################################################
################### SourMash Analysis ###################
#########################################################

# https://sourmash.readthedocs.io/en/latest/command-line.html#sourmash-compare-compare-many-signatures
# https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html#parameter-strings

rule Sourmash_Sketch_Genome_Scaled1To1:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
    output:
        Sourmash_Sketch_Asm_SIG = output_Dir + "/AsmAnalysis/SourMash/{sampleID}.LR.Asm.SourMash.Scaled1.sig",
    conda:
        "CondaEnvs/sourmash_v4_4_3.yml"
    threads: 1
    shell:
        "sourmash sketch dna {input.i_Assembly_FA} -p scaled=1 -o {output.Sourmash_Sketch_Asm_SIG} "

rule Sourmash_Compare_Scaled1To1:
    input:
        expand(output_Dir + "/AsmAnalysis/SourMash/{sampleID}.LR.Asm.SourMash.Scaled1.sig", sampleID = input_All_SampleIDs),       
    output:
        Sourmash_Compare_TXT = output_Dir + "/SourMash/CompareAllAsm_Scaled1/SourMash.Compare.Scaled1.out",
        Sourmash_Compare_Plot_Matrix_PDF = output_Dir + "/SourMash/CompareAllAsm_Scaled1/SourMash.Compare.Scaled1.out.matrix.pdf",
        Sourmash_Compare_Plot_Dendro_PDF = output_Dir + "/SourMash/CompareAllAsm_Scaled1/SourMash.Compare.Scaled1.out.dendro.pdf",
    conda:
        "CondaEnvs/sourmash_v4_4_3.yml"
    threads: 1
    params:
        Compare_OutDir = output_Dir + "/SourMash/CompareAllAsm_Scaled1/"
    shell:
        "sourmash compare {input} -o {output.Sourmash_Compare_TXT} \n"
        "sourmash plot --pdf --labels {output.Sourmash_Compare_TXT} --output-dir {params.Compare_OutDir} "

#########################################################















