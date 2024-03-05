

### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]

H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]
H37rv_GBK_PATH = config["H37rv_GBK_PATH"]


Kraken2_DB_PATH = config["Kraken2_DB_PATH"]


# Define PATH of main OUTPUT directory
output_Dir = config["output_dir"]


# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')



# Create a python list of Sample IDs
input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )
input_SampleIDs_SRA_acc = list( input_DataInfo_DF["Short Reads Accession"].values )


SampleID_to_SRARunAcc_Dict = dict(input_DataInfo_DF[['SampleID', 'Short Reads Accession']].values)

print("List of input sampleIDs:", len(input_All_SampleIDs), input_All_SampleIDs)
print("List of input SRA Run Accs:", len(input_SampleIDs_SRA_acc), input_SampleIDs_SRA_acc)


rule all:
    input:
        expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gbff", sampleID_WiIll = input_All_SampleIDs),

        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",


        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
        output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",

        output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",

        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",



# Step 4: Annotate SR de novo asm w/ Bakta

rule Bakta_Anno_SRAsm:
    input:
        i_Assembly_FA = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.fasta",
    output:
        Asm_Bakta_Anno_FAA = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.faa",
        Asm_Bakta_Anno_FNA = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.fna",
        Asm_Bakta_Anno_GBFF = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gbff",
        Asm_Bakta_Anno_GFF = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3",
    conda:
        "CondaEnvs/Bakta_1_6_1.nobuilds.yml" 
    threads: 8
    params:
        Bakta_OutputDir_PATH = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/",
        Bakta_DB_Dir = "/n/data1/hms/dbmi/farhat/mm774/References/Bakta_DB_v4/db",
        Bakta_OutPrefix = "{sampleID_WiIll}.Bakta"
    shell:
        "bakta --db {params.Bakta_DB_Dir} --verbose --output {params.Bakta_OutputDir_PATH} "
        " --prefix {params.Bakta_OutPrefix} --locus-tag {wildcards.sampleID_WiIll}  "
        " --threads {threads} {input.i_Assembly_FA} "


####################################################################################





##################################################
############ Pangenome Analysis ##################
##################################################

rule Panaroo_StdParams_Strict_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes"


rule Panaroo_StdParams_Moderate_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes"


rule Panaroo_StdParams_Sensitive_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes"




##### Run Panaroo but with the '--merge_paralogs' option. 

rule Panaroo_StdParams_Strict_MergeParalogs_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Moderate_MergeParalogs_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Sensitive_MergeParalogs_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/summary_statistics.txt",
       gene_data_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_data.csv",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_Strict_MergeParalogs_AvA_KmerComparison_SRAsm:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "

rule Panaroo_Moderate_MergeParalogs_AvA_KmerComparison_SRAsm:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp,
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "

rule Panaroo_Sensitive_MergeParalogs_AvA_KmerComparison_SRAsm:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "







rule Roary_DefaultParams_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        'TMPDIR={params.Roary_TmpDir}  \n'
        'mkdir -p {params.Roary_TmpDir}  \n'
        'echo $TMPDIR \n'
        "roary -e -n -v -r -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e --mafft


rule Roary_NoSplitParalogs_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        'TMPDIR={params.Roary_TmpDir}  \n'
        'mkdir -p {params.Roary_TmpDir}  \n'
        'echo $TMPDIR \n'
        "roary -e -n -v -r -s -p {threads} -f {params.Roary_OutputDir} {input} " # -s means 'dont split paralogs'
        # -e --mafft

rule Roary_NoSplitParalogs_I90_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        'TMPDIR={params.Roary_TmpDir}  \n'
        'mkdir -p {params.Roary_TmpDir}  \n'
        'echo $TMPDIR \n'
        "roary -e -n -v -r -i 90 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -s means 'dont split paralogs'
        # -i 90 sets the minimumblastp identity to 90%
        # -e --mafft
 

rule Roary_NoSplitParalogs_I80_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/GenomeAnnotation/{sampleID_WiIll}_Asm_Bakta/{sampleID_WiIll}.Bakta.gff3", sampleID_WiIll = input_All_SampleIDs),           
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/",
        Roary_TmpDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/tmp",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        'TMPDIR={params.Roary_TmpDir}  \n'
        'mkdir -p {params.Roary_TmpDir}  \n'
        'echo $TMPDIR \n'
        "roary -e -n -v -r -i 80 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -s means 'dont split paralogs'
        # -i 80 sets the minimumblastp identity to 80%
        # -e --mafft
 





rule Roary_DefaultParams_AvA_KmerComparison_SRAsm:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "

rule Roary_NoSplitParalogs_AvA_KmerComparison_SRAsm:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "


rule Roary_NoSplitParalogs_I90_AvA_KmerComparison_SRAsm:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "



rule Roary_NoSplitParalogs_I80_AvA_KmerComparison_SRAsm:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time pgqc ava -k {params.k} -i {input} -o {output} "





###### Run GeneSeqCheck step for all Panaroo & Roary analyses ######

# i_Assembly_FA = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.fasta",




rule create_SampleToAsmFA_TSV_SRAsm:
    output:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
    run:
        
        # Template file path with a placeholder for SampleID
        base_asm_path = output_Dir + '/SR_DataProcessing/{}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.fasta'

        input_DataInfo_DF['ShortRead_Genome_ASM_PATH'] = input_DataInfo_DF['SampleID'].apply(lambda x: base_asm_path.format(x))

        # Output as tsv
        input_DataInfo_SRAsm_DF = input_DataInfo_DF[['SampleID', 'ShortRead_Genome_ASM_PATH']].copy()
        input_DataInfo_SRAsm_DF.columns = ['SampleID', 'Genome_ASM_PATH']
        input_DataInfo_SRAsm_DF.to_csv(output.asmfa_tsv, sep="\t", index = False)


rule Panaroo_Strict_MergeParalogs_PGQC_AsmGeneSeqCheck_SRAsm:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Panaroo_Moderate_MergeParalogs_PGQC_AsmGeneSeqCheck_SRAsm:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Panaroo_Sensitive_MergeParalogs_PGQC_AsmGeneSeqCheck_SRAsm:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Roary_DefaultParams_PGQC_AsmGeneSeqCheck_SRAsm:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",  
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Roary_NoSplitParalogs_PGQC_AsmGeneSeqCheck_SRAsm:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Roary_NoSplitParalogs_I90_PGQC_AsmGeneSeqCheck_SRAsm:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Roary_NoSplitParalogs_I80_PGQC_AsmGeneSeqCheck_SRAsm: 
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time pgqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


#####################################################################















































