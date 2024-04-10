# 5.Ecoli.Generate.SRAsms.smk
### Snakemake - Run pipeline for short-read de novo assembly of 50 E. coli isolates (From Shaw et. al. 2021)
### Maximillian Marin (mgmarin@g.harvard.edu)

### Import Statements ###
import pandas as pd


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


############################# SR WGS Assembly & Analysis #############################

# Step 1: Download Paired-end SR WGS associated w/ 
rule download_FQ_FromSRA_RunID:
    output: 
         fq1_unzipped = temp(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq"),
         fq2_unzipped = temp(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq"),
    params:
        target_DownloadDir = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/",
        input_SRA_Acc = lambda wildcards: SampleID_to_SRARunAcc_Dict[wildcards.sampleID_WiIll],
    conda:
        "CondaEnvs/sratools_2_10_7_Conda.yml"
    shell:
        "fastq-dump --split-files --outdir {params.target_DownloadDir} {params.input_SRA_Acc}\n"
        "mv {params.target_DownloadDir}/{params.input_SRA_Acc}_1.fastq {output.fq1_unzipped} \n"
        "mv {params.target_DownloadDir}/{params.input_SRA_Acc}_2.fastq {output.fq2_unzipped} \n"


rule GZIP_Illumina_FQs:
    input:
        fq1_unzipped = (output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq"),
        fq2_unzipped = (output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq"),
    output:
        fq1_gz = temp(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq.gz"),
        fq2_gz = temp(output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq.gz"),
    shell:
        "gzip -c < {input.fq1_unzipped} > {output.fq1_gz} \n"
        "gzip -c < {input.fq2_unzipped} > {output.fq2_gz}"


# Step 2: Trim reads w/ FASTP

rule fastp_Illumina_PE_Trim:
    input:
        fq1 = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_1.fastq.gz",
        fq2 = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs/{sampleID_WiIll}_2.fastq.gz",
    output:
        fq1_trimmed = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiIll}_1.fastp.trim.fastq.gz",
        fq2_trimmed = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiIll}_2.fastp.trim.fastq.gz",
        html = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiIll}.fastp.html",
        json = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiIll}.fastp.json",
    conda: "CondaEnvs/fastp_v0_23_3.yml"
    params:
        min_length = "50"
    shell:
        "fastp -i {input.fq1} -I {input.fq2} -o {output.fq1_trimmed} -O {output.fq2_trimmed} -h {output.html} -j {output.json} --length_required {params.min_length}"



# Step 3: Run SPAdes assembler

# Adding SPADes assembly with UNIcycler to analysis

### A) Assembly with SPAdes through Unicycler
rule unicycler_SPAdes_Assemble_IlluminaWGS:
    input:
        fq1_trimmed = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiIll}_1.fastp.trim.fastq.gz",
        fq2_trimmed = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiIll}_2.fastp.trim.fastq.gz",
        #DnaA_Seq_fa = "references/DnaA_MTb_H37Rv_dna.fasta"   #H37rv_DnaA_FA_PATH,
    output:
        assembly_GFA = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.gfa",
        assembly_fa = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/assembly.fasta",
    conda: "CondaEnvs/unicycler_4_8.nobuilds.yml"
    threads: 8
    params:
        Unicycler_OutputDir_PATH = output_Dir + "/SR_DataProcessing/{sampleID_WiIll}/IlluminaWGS/Unicycler_SPAdesAssembly/"
    shell:
        "unicycler  -t {threads} " #--start_genes {input.DnaA_Seq_fa} "
        " -1 {input.fq1_trimmed} -2 {input.fq2_trimmed} "
        " -o {params.Unicycler_OutputDir_PATH} "
        
        # " ‑‑mode conservative " # This version of unicycler doesn't support the --mode arguement




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










