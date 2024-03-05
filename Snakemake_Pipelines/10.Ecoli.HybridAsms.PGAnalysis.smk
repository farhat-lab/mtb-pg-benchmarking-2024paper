
# PMP_SM - Preprocessing for merging PacBio FQs (Any # of FQs to single merged FQ for one isolate)


# This is a test run of the merging approach for FQs


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


rule all:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.fna", sampleID = input_All_SampleIDs),

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/summary_statistics.txt",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/summary_statistics.txt",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/summary_statistics.txt",

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   

        
        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",


        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   

        output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv", 

        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",

        output_Dir + "/Phylogenies/iqtree/Panaroo.Strict.MergeParalogs.MidRoot.WiNodeNames.newick",
        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/Panstripe_PhyloAnalysis/Panstripe.PanstripeSummary.tsv",


        output_Dir + "/FastANI/FastANI_LRAsm/FastANI.AllVsAll.LRAsm.txt",

        output_Dir + "/SourMash/CompareAllAsm_Default/SourMash.Compare.Default.out",
        output_Dir + "/SourMash/CompareAllAsm_Scaled1/SourMash.Compare.Scaled1.out",




##################################################
############ Genome Assembly QC ##################
##################################################


# Run BUSCO for general QC and completeness evaluation

# ......


rule busco_Ecoli_LRAsm:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
    output:
        ShortSum_TXT = output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_LRAsm/{sampleID}/short_summary.specific.enterobacterales_odb10.{sampleID}.txt",
    conda:
        "CondaEnvs/Busco_v5_4_4.nobuilds.yml"
    params:
        mode = "genome",
        name = '{sampleID}',
        out_dir = output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_LRAsm/",
        tmp_dir = output_Dir + "/Busco_Download_Tmp/",
        lineage="enterobacterales_odb10",
    shell:
        #"mkdir {params.tmp_dir} \n"
        "mkdir -p {params.out_dir} \n"
        "cd {params.out_dir} \n"
        "busco -i {input.i_Assembly_FA} -o {wildcards.sampleID} --out_path {params.out_dir} --mode {params.mode} "
        " --lineage_dataset {params.lineage} --download_path {params.tmp_dir} --force " #--offline "






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


# Relevant link: https://github.com/gtonkinhill/panaroo/issues/68
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





##################################################
############ Pangenome Analysis ##################
##################################################

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


rule Panaroo_StdParams_Strict_MergeParalogs_AllIsolates_WiCoreAln:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/pan_genome_reference.fa",
       iq_tree = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/core_tree.treefile"
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 8
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo --alignment core -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes --merge_paralogs \n"
        "\n"
        "cd {params.Panaroo_OutputDir} \n"
        "iqtree -s core_gene_alignment.aln -pre core_tree -nt 1 -fast -m GTR "


rule IQtree_GTR_Panaroo_S_MP_CoreAln:
    input:
       CoreAln_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/core_gene_alignment.aln"
    output:
        TREE = output_Dir + "/Phylogenies/iqtree/Panaroo.Strict.MergeParalogs.treefile",
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree/Panaroo.Strict.MergeParalogs",
        bb = 10000
    shell:
        "time iqtree -m GTR -nt 1 -bb {params.bb} -s {input} -pre {params.iqtree_outprefix} " # -m # -redo


rule Update_IQtree_Newick_10AmbPLCMask_Filt:
    input:
        TREE = output_Dir + "/Phylogenies/iqtree/Panaroo.Strict.MergeParalogs.treefile",
    output:
        TREE_WiNodeNames = output_Dir + "/Phylogenies/iqtree/Panaroo.Strict.MergeParalogs.MidRoot.WiNodeNames.newick",
    shell:
        "Scripts/IQTree.UpdateNodeNames.V1.py -i {input} -o {output}"


##### Run Panstripe Phylo Analysis w/ GAIN + LOSS measure #####

rule RunPanstripe_Panaroo_S_MP_WiCoreAln:
    input:
       PA_Rtab = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/gene_presence_absence.Rtab",
       in_tree = output_Dir + "/Phylogenies/iqtree/Panaroo.Strict.MergeParalogs.MidRoot.WiNodeNames.newick",
    output:
       PA_NoBakta_Rtab = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/gene_presence_absence.Renamed.Rtab",
       Panstripe_Summ_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/Panstripe_PhyloAnalysis/Panstripe.PanstripeSummary.tsv"
    threads: 1
    params:
        Panstripe_OutputDir = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates_WiCoreAln/Panstripe_PhyloAnalysis",
    shell:
        'rm -r {params.Panstripe_OutputDir} \n'
        "~/conda3/bin/conda deactivate \n"
        "module load gcc/6.2.0 R/4.1.1 \n"
        'export PATH="/n/data1/hms/dbmi/farhat/mm774/Snakemake_Pipelines/Mtb-WGA-SMK/Scripts/Panstripe_Scripts:$PATH" \n'
        ""
        " sed 's/\.Bakta//g' {input.PA_Rtab} > {output.PA_NoBakta_Rtab} \n"
        "RunPhyloAnalysis.WiPanstripe.V1.R -p {input.in_tree} -r {output.PA_NoBakta_Rtab} -o {params.Panstripe_OutputDir} -a Panstripe"
















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
       #SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/summary_statistics.txt",
       #gene_data_CSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary -e -n -r -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e --mafft


rule Roary_NoSplitParalogs_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       #SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/summary_statistics.txt",
       #gene_data_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary -e -n -r -s -p {threads} -f {params.Roary_OutputDir} {input} " # -s means 'dont split paralogs'
        # -e --mafft



rule Roary_NoSplitParalogs_I90_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       #SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/summary_statistics.txt",
       #gene_data_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary -e -n -r -i 90 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -s means 'dont split paralogs'
        # -i 90 sets the minimumblastp identity to 90%
        # -e --mafft



rule Roary_NoSplitParalogs_I80_AllIsolates:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
    output:
       #SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/summary_statistics.txt",
       #gene_data_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary -e -n -r -i 80 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -s means 'dont split paralogs'
        # -i 80 sets the minimumblastp identity to 80%
        # -e --mafft
 









rule Panaroo_Strict_MergeParalogs_AvA_KmerComparison:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -i {input} -o {output} "
        #"time Scripts/KmerComp.AllvsAll.V3.py -k {params.k} -i {input} -o {output}  "

rule Panaroo_Moderate_MergeParalogs_AvA_KmerComparison:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp,
    shell:
        "time panqc ava -i {input} -o {output} "
        #"time Scripts/KmerComp.AllvsAll.V3.py -k {params.k} -i {input} -o {output} "

rule Panaroo_Sensitive_MergeParalogs_AvA_KmerComparison:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp 
    shell:
        "time panqc ava -i {input} -o {output} "
        #"time Scripts/KmerComp.AllvsAll.V3.py -k {params.k} -i {input} -o {output} "




rule Roary_DefaultParams_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -i {input} -o {output} "

rule Roary_NoSplitParalogs_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -i {input} -o {output} "


rule Roary_NoSplitParalogs_I90_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -i {input} -o {output} "



rule Roary_NoSplitParalogs_I80_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "



###### Run GeneSeqCheck step for all Panaroo & Roary analyses ######

rule create_SampleToAsmFA_TSV:
    output:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.FA.PATHs.tsv",
    run:
        # Output as tsv
        input_DataInfo_DF[['SampleID', 'Genome_ASM_PATH']].to_csv(output.asmfa_tsv, sep="\t", index = False)


rule Panaroo_Strict_MergeParalogs_PGQC_AsmGeneSeqCheck:
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

rule Panaroo_Moderate_MergeParalogs_PGQC_AsmGeneSeqCheck:
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



rule Panaroo_Sensitive_MergeParalogs_PGQC_AsmGeneSeqCheck:
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




rule Roary_DefaultParams_PGQC_AsmGeneSeqCheck:
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

rule Roary_NoSplitParalogs_PGQC_AsmGeneSeqCheck:
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


rule Roary_NoSplitParalogs_I90_PGQC_AsmGeneSeqCheck:
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



rule Roary_NoSplitParalogs_I80_PGQC_AsmGeneSeqCheck: 
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







#########################################################
############ LR Asm FastA ANI Analysis ##################
#########################################################


rule Save_LR_Asm_PATH_To_TXT: 
    input:
        Asm_Bakta_Anno_FNA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_Asm_Bakta/{sampleID}.Bakta.fna",
    output:
        output_Dir + "/AsmAnalysis/{sampleID}/FastANI/{sampleID}.PathToFASTA.LRAsm.txt"
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


# https://sourmash.readthedocs.io/en/latest/sourmash-sketch.html#parameter-strings

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



