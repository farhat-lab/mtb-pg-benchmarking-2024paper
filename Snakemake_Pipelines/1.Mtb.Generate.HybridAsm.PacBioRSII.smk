# 1.Mtb.Generate.HybridAsm.PacBioRSII.smk
### Snakemake - Run pipeline for PacBio subread (RSII) long read de novo assembly (+ short-read polishing of base-level errors)


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



input_DataInfo_DF_With_PacBio_WGS = input_DataInfo_DF[  input_DataInfo_DF["PacBio_FQ_PATH"] != "None" ]

input_DataInfo_DF_With_Illumina_WGS = input_DataInfo_DF[ input_DataInfo_DF["Illumina_PE_FQs_PATH"] != "None" ]

input_DataInfo_DF_With_Both_Illumina_And_PacBio_WGS = input_DataInfo_DF[ (input_DataInfo_DF["Illumina_PE_FQs_PATH"] != "None") & (input_DataInfo_DF["PacBio_FQ_PATH"] != "None") ]




# Create a python list of Sample IDs

input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )

input_SampleIDs_With_PB_And_Ill = list( input_DataInfo_DF_With_Both_Illumina_And_PacBio_WGS["SampleID"].values )


SampleIDTo_Illumina_PE_FQ1_Dict = {}
SampleIDTo_Illumina_PE_FQ2_Dict = {}

SampleID_To_PacBio_FQs_Dict = {}


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

        
    if SampleID_i in input_SampleIDs_WithPacBio: # If PacBio WGS is provided

        PacBio_FQs_PATH = row["PacBio_FQ_PATH"]
        PacBio_FastQ_Files_List = PacBio_FQs_PATH.split(";")

        SampleID_To_PacBio_FQs_Dict[SampleID_i] = PacBio_FastQ_Files_List



print("SampleIDs with PacBio WGS:", len(input_SampleIDs_WithPacBio), input_SampleIDs_WithPacBio)
print("SampleIDs with Illumina WGS:", len(input_SampleIDs_WithIllumina), input_SampleIDs_WithIllumina)






rule all:
    input:
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz", sampleID_WiPacBio = input_SampleIDs_WithPacBio),  
        expand(output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly.fasta", sampleID_WiPacBio = input_SampleIDs_WithPacBio),
        expand(output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta.fai", sampleID_WiPacBio = input_SampleIDs_WithPacBio),



rule merge_Wi_CAT_All_PacBio_FQ_GZs:
    input:
        lambda wildcards: expand("{i_PacBio_FQ_PATH}", i_PacBio_FQ_PATH = SampleID_To_PacBio_FQs_Dict[wildcards.sampleID_WiPacBio])
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
    threads: 1
    shell:
        "cat {input} > {output}"



###################################################
##### Flye (create assembly & polish) ######
###################################################


rule flye_Assemble: # Flye v2.6 w/ asmCov = 200
    input:
        pb_subreads_fq = output_Dir + "/{sampleID_WiPacBio}/pacbio/pacbio_reads/{sampleID_WiPacBio}.merged.subreads.fastq.gz",
    output:
        assembly_fa = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly.fasta",
        assembly_info_txt = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly_info.txt"
    conda:
       "envs/PacBio_Software_py27_Conda.yml"
    threads: 10
    params:
        Flye_OutputDir_PATH = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/"
    shell:
        "flye --pacbio-raw {input.pb_subreads_fq} --out-dir {params.Flye_OutputDir_PATH} "
        "--genome-size 5m --threads {threads} --asm-coverage 200 --iterations 3"


###################################################################################
######### CIRCLATOR for setting start at DnaA (Assuming Circular genome) ##########
###################################################################################

rule circlator_FixStart_DnaA:
    input:
        #assembly_not_circularcontigs_txt = output_Dir + "{sampleID}/pacbio/Flye_Assembly/NoncircularContigs_List.txt",
        DnaA_Seq_fa = H37rv_DnaA_FA_PATH,
        flye_assembly_fa = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/assembly.fasta"
    output:
        flye_assembly_FixStart_assembly = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta"
    conda:
        "envs/circulator_151py37_Conda.yml"
    threads: 1
    shell:
        "circlator fixstart --genes_fa {input.DnaA_Seq_fa} {input.flye_assembly_fa} {output_Dir}/{wildcards.sampleID_WiPacBio}/pacbio/Flye_Assembly/{wildcards.sampleID_WiPacBio}.flyeassembly.fixstart"





rule samtools_faidx_FlyeAssembly:
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta"
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta.fai"
    conda:
        "envs/PacBio_Software_py27_Conda.yml"
    threads: 1
    shell: "samtools faidx {input}"




rule CP_PacBio_FlyeAssembly_To_I3_Dir: # This is a TEMP fix to just skip the quiver polishing steps
    input:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly/{sampleID_WiPacBio}.flyeassembly.fixstart.fasta"
    output:
        output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.fixstart.fasta",
    shell: 
        "cp {input} {output}"




rule filterByLength_100kbContigs_FlyeAssembly_I3:
    input:
        PacBio_Flye_Assembly_fa = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.fixstart.fasta",
    output:
        PacBio_Flye_Assembly_Renamed_100KbContigs_FA = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta",
        PacBio_Flye_Assembly_Renamed_100KbContigs_FAI = output_Dir + "/{sampleID_WiPacBio}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiPacBio}.flyeassembly.I3.Renamed.100Kb.fasta.fai"
    threads: 1
    conda:
        "envs/bioinfo_util_env_V1.yml" # "envs/bioinfo_util_env.yml"
    shell:
        " bioawk -c fastx '{{ print \">{wildcards.sampleID_WiPacBio}_\"$name \"\\n\" $seq }}' {input.PacBio_Flye_Assembly_fa} "
        " | "
        " bioawk -c fastx '{{ if(length($seq) > 100000) {{ print \">\"$name; print $seq }}}}' > {output.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} \n"
        " samtools faidx {output.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} "







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




##### Combined Illumina + PacBio Analysis Steps #####



rule bwa_idx_ref_100kbContigs_FlyeAssembly_I3:
    input:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta.bwt"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "bwa index {input}"


rule bwa_map_IllPE_AlignTo_I3_Assembly:
    input:
        PacBio_Flye_Assembly_I3_FA = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta",
        PacBio_Flye_Assembly_I3_FA_BWT_IDX = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta.bwt",
        fq1_trimmed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_Wi_Ill_And_PB}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_Wi_Ill_And_PB}/IlluminaWGS/FASTQs_Trimmomatic_Trimming/{sampleID_Wi_Ill_And_PB}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.sam"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_Wi_Ill_And_PB}\tSM:{sampleID_Wi_Ill_And_PB}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.PacBio_Flye_Assembly_I3_FA} "
        "{input.fq1_trimmed} {input.fq2_trimmed} > {output}"


rule samtools_ViewAndSort_IllPE_AlignTo_I3_Assembly:
    input:
        output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.sam"
    output:
        bam = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.bam",
        bai = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.bam.bai",
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} | samtools sort - > {output.bam} \n"
        "samtools index {output.bam}"



#####################################
#### PICARD (remove duplicates) #####
#####################################

rule picard_RemoveDup_IllPE_AlignTo_I3_Assembly:
    input:
        IllPE_BwaMEM_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.bam",
    output:
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam",
        IllPE_BwaMEM_Duprem_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam.bai",
        IllPE_BwaMEM_Duprem_METRICS = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam.metrics",
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "picard -Xmx5g MarkDuplicates I={input.IllPE_BwaMEM_BAM} O={output.IllPE_BwaMEM_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.IllPE_BwaMEM_Duprem_METRICS} ASSUME_SORT_ORDER=coordinate"
        " \n"
        "samtools index {output.IllPE_BwaMEM_Duprem_BAM}"




rule pilon_IllPE_Polishing_I3_Assembly:
    input:
        PacBio_Flye_Assembly_I3_FA = output_Dir + "/{sampleID_Wi_Ill_And_PB}/pacbio/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_Wi_Ill_And_PB}.flyeassembly.I3.Renamed.100Kb.fasta",
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam",
        IllPE_BwaMEM_Duprem_BAI = output_Dir + "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/IlluminaPE_AlignedTo_Q3Assembly/{sampleID_Wi_Ill_And_PB}.IllPE.AlnTo.I3Assembly.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.vcf",
        pilon_I3_Polished_FA = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
        pilon_I3_PP_ChangesFile = output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.changes"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx10g --fix snps,indels --genome {input.PacBio_Flye_Assembly_I3_FA} --bam {input.IllPE_BwaMEM_Duprem_BAM} --output {wildcards.sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished"
        " --outdir {output_Dir}/{wildcards.sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/ --variant --changes --tracks \n"
        "samtools faidx {output.pilon_I3_Polished_FA}"



rule samtools_faidx_PilonPolished_I3Assembly:
    input:
        output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta",
    output:
        output_Dir+ "/{sampleID_Wi_Ill_And_PB}/FlyeAssembly_I3_IlluminaPolishing/pilon_IllPE_Polishing_I3_Assembly_ChangeSNPsINDELsOnly/{sampleID_Wi_Ill_And_PB}.Flye.I3Assembly.PilonPolished.fasta.fai"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    threads: 1
    shell:
        "samtools faidx {input}"




















