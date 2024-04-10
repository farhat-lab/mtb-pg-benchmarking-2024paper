# 4.Mtb.Generate.SRAsms.smk
### Snakemake - Run pipeline for short-read de novo assembly for Mtb isolates
### Maximillian Marin (mgmarin@g.harvard.edu)

### Import Statements ###
import pandas as pd

### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]

H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]

Kraken2_DB_PATH = config["Kraken2_DB_PATH"]


# Define PATH of main OUTPUT directory
output_Dir = config["output_dir"]



# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')
input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )


input_DataInfo_DF_With_Illumina_WGS = input_DataInfo_DF[ input_DataInfo_DF["Illumina_PE_FQs_PATH"] != "None" ]
input_SampleIDs_WithIllumina = list( input_DataInfo_DF_With_Illumina_WGS["SampleID"].values )

SampleIDTo_Illumina_PE_FQ1_Dict = {}
SampleIDTo_Illumina_PE_FQ2_Dict = {}

# Iterate over each sample's row and define path to input FQ files for each sequencing technology

for idx, row in input_DataInfo_DF.iterrows():
    
    SampleID_i = row["SampleID"]
    
    if SampleID_i in input_SampleIDs_WithIllumina: # If Illumina WGS is provided

        Illumina_FQs_PATH = row["Illumina_PE_FQs_PATH"]
        Ill_FastQ_Files_List = Illumina_FQs_PATH.split(";")
        

        if len(Ill_FastQ_Files_List) == 2:

            FQ_1_PATH, FQ_2_PATH = Ill_FastQ_Files_List
            SampleIDTo_Illumina_PE_FQ1_Dict[SampleID_i] = FQ_1_PATH
            SampleIDTo_Illumina_PE_FQ2_Dict[SampleID_i] = FQ_2_PATH


print("SampleIDs with Illumina WGS:", len(input_SampleIDs_WithIllumina), input_SampleIDs_WithIllumina)







rule all:
    input:
        expand(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.fastq", sampleID_WiONT = input_SampleIDs_WithONT),                     
        expand(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/assembly.fasta", sampleID_WiONT = input_SampleIDs_WithONT),  
        expand(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/{sampleID_WiONT}.SRAsm.Kraken2.MTBConly.fasta", sampleID_WiONT = input_SampleIDs_WithONT), 



# adapter list from: https://github.com/stephenturner/adapters/blob/master/adapters_combined_256_unique.fasta
# Also adapter list combined w/ ftp://ftp.ncbi.nlm.nih.gov/pub/kitts/adaptors_for_screening_proks.fa
# Adapter list path: references/CustomTrimmomatic_IlluminaWGS_AdapterList.WiProkAdaptersNCBI.fasta

rule trimmomatic_Illumina_PE_Trimming_V2:
    input:
        r1 = lambda wildcards: SampleIDTo_Illumina_PE_FQ1_Dict[wildcards.sampleID_WiONT],
        r2 = lambda wildcards: SampleIDTo_Illumina_PE_FQ2_Dict[wildcards.sampleID_WiONT],
    output:
        r1 = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.fastq",
        r2 = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_2_trimmed.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.unpaired.fastq",
        r2_unpaired = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_2_trimmed.unpaired.fastq",
    log:
        output_Dir + "/logs/trimmomatic/{sampleID_WiONT}.v2.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:'references/CustomTrimmomatic_IlluminaWGS_AdapterList.WiProkAdaptersNCBI.fasta':2:30:10:2:true SLIDINGWINDOW:4:20 MINLEN:75"],
        # optional parameters
        # extra=" "
    threads: 8
    wrapper:
        "0.38.0/bio/trimmomatic/pe"



#### run Kraken on PE short reads ####

rule Kraken2_Illumina_PE:
    input:
        Kraken2_DB = Kraken2_DB_PATH,
        #fq1_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiONT}_1.fastp.trim.fastq.gz",
        #fq2_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiONT}_2.fastp.trim.fastq.gz",
        fq1_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_2_trimmed.fastq",
    output:
        ReadClassifications = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Kraken2.Output.Reads.txt",
        Kraken2_Report = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Kraken2.Output.Report.txt"
    conda:
        "envs/kraken2_2_0_8_Conda.yml"
    threads: 1
    shell:
        " kraken2 --use-names --threads {threads} --db {input.Kraken2_DB} --output {output.ReadClassifications} "
        " --report {output.Kraken2_Report} --paired {input.fq1_trimmed} {input.fq2_trimmed}  "




# https://github.com/jenniferlu717/KrakenTools#extract_kraken_readspy 
rule Filt_Reads_MTBCandAbove_Only:
    input:
        #fq1_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiONT}_1.fastp.trim.fastq.gz",
        #fq2_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Fastp_Trimming/{sampleID_WiONT}_2.fastp.trim.fastq.gz",
        fq1_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_2_trimmed.fastq",
        ReadClassifications = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Kraken2.Output.Reads.txt",
        Kraken2_Report = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Kraken2.Output.Report.txt"
    output:
        R1_KrakFilt_FQ = temp(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Trimmed.MTBC.Filtered.R1.fastq"),
        R2_KrakFilt_FQ = temp(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Trimmed.MTBC.Filtered.R2.fastq"),
        R1_KrakFilt_FQ_GZ = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Trimmed.MTBC.Filtered.R1.fastq.gz",
        R2_KrakFilt_FQ_GZ = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Trimmed.MTBC.Filtered.R2.fastq.gz",
    conda:
        "envs/krakentools_1_2.yml"
    shell:
        "extract_kraken_reads.py -k {input.ReadClassifications} "
        " --report {input.Kraken2_Report} "
        " -s {input.fq1_trimmed} -s2 {input.fq2_trimmed}"
        " -o {output.R1_KrakFilt_FQ} -o2 {output.R2_KrakFilt_FQ} "
        " --fastq-output "
        " -t 77643 --include-children --include-parents > /dev/null \n" # reads classified as MTBC (ABOVE & BELOW) will be extracted
        "gzip -c {output.R1_KrakFilt_FQ} > {output.R1_KrakFilt_FQ_GZ} \n"
        "gzip -c {output.R2_KrakFilt_FQ} > {output.R2_KrakFilt_FQ_GZ} \n"




### Let's run SPAdes assembler

# Adding SPADes assembly with UNIcycler to analysis

### A) Assembly with SPAdes through Unicycler v 4.8
rule unicycler_v048_SPAdes_SR_Asm_KrakFiltReads:
    input:
        R1_KrakFilt_FQ = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Trimmed.MTBC.Filtered.R1.fastq.gz",
        R2_KrakFilt_FQ = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2/{sampleID_WiONT}.Trimmed.MTBC.Filtered.R2.fastq.gz",
        DnaA_Seq_fa = "references/DnaA_MTb_H37Rv_dna.fasta"   #H37rv_DnaA_FA_PATH,
    output:
        assembly_GFA = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/assembly.gfa",
        assembly_fa = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/assembly.fasta",
    conda:
        "envs/unicycler_4_8.yml"
    threads: 8
    params:
        Unicycler_OutputDir_PATH = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/"
    shell:
        "rm -r {params.Unicycler_OutputDir_PATH}/ \n"
        "unicycler -t {threads}  --start_genes {input.DnaA_Seq_fa} "
        " -1 {input.R1_KrakFilt_FQ} -2 {input.R2_KrakFilt_FQ} "
        " -o {params.Unicycler_OutputDir_PATH} "



rule Kraken2_SRAsm:
    input:
        Kraken2_DB = Kraken2_DB_PATH,
        i_Assembly_FA = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/assembly.fasta",
    output:
        ReadClassifications = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2_SRAsm/{sampleID_WiONT}.Kraken2.Reads.txt",
        Kraken2_Report      = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2_SRAsm/{sampleID_WiONT}.Kraken2.Report.txt"
    conda:
        "envs/kraken2_2_0_8_Conda.yml"
    threads: 1
    shell:
        " kraken2 --use-names --threads {threads} --db {input.Kraken2_DB} --output {output.ReadClassifications} "
        " --report {output.Kraken2_Report} {input.i_Assembly_FA} "




rule Filt_SRAsm_Contigs_MTBC_Only:
    input:
        i_Assembly_FA = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/assembly.fasta",
        ReadClassifications = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2_SRAsm/{sampleID_WiONT}.Kraken2.Reads.txt",
        Kraken2_Report      = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Kraken2_SRAsm/{sampleID_WiONT}.Kraken2.Report.txt"
    output:
        i_Assembly_KrakFilt_FA = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/{sampleID_WiONT}.SRAsm.Kraken2.MTBConly.fasta",
    conda:
        "envs/krakentools_1_2.yml"
    shell:
        "extract_kraken_reads.py -k {input.ReadClassifications} --report {input.Kraken2_Report} -s {input.i_Assembly_FA} -o {output} -t 77643 --include-children"  # reads classified as MTBC (& BELOW) will be extracted













