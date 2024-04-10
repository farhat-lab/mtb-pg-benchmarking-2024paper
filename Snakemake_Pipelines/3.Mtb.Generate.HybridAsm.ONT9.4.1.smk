# 3.Mtb.Generate.HybridAsm.ONT9.4.1.smk
### Snakemake - Run pipeline for ONT (v9.4.1) long read de novo assembly (+ short-read polishing of base-level errors)
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


input_DataInfo_DF_With_ONT_WGS = input_DataInfo_DF[  input_DataInfo_DF["ONT_FQ_PATH"] != "None" ]

SampleID_To_ONT_FQ_Dict = dict(input_DataInfo_DF_With_ONT_WGS[["SampleID", "ONT_FQ_PATH"]].values)
input_SampleIDs_WithONT = list( input_DataInfo_DF_With_ONT_WGS["SampleID"].values )



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



print("SampleIDs with ONT WGS:", len(input_SampleIDs_WithONT), input_SampleIDs_WithONT)
print("SampleIDs with Illumina WGS:", len(input_SampleIDs_WithIllumina), input_SampleIDs_WithIllumina)


rule all:
    input:
        expand(output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly/assembly.fasta", sampleID_WiONT = input_SampleIDs_WithONT),
        expand(output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiONT}.flyeassembly.I3.Renamed.100Kb.fasta", sampleID_WiONT = input_SampleIDs_WithONT),
        expand(output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta", sampleID_WiONT = input_SampleIDs_WithONT),
        
        expand(output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.fasta", sampleID_WiONT = input_SampleIDs_WithONT),  
        expand(output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.fixstart.fasta", sampleID_WiONT = input_SampleIDs_WithONT), 

        expand(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/assembly.fasta", sampleID_WiONT = input_SampleIDs_WithONT),  
        
        expand(output_Dir + "/{sampleID_WiONT}/IlluminaWGS/Unicycler_SPAdesAssembly_v0_4_8/{sampleID_WiONT}.SRAsm.Kraken2.MTBConly.fasta", sampleID_WiONT = input_SampleIDs_WithONT), 






rule flye_Assemble_ONTraw: # Flye v2.6 w/ asmCov = 200
    input:
        ont_reads_fq = lambda wildcards: SampleID_To_ONT_FQ_Dict[wildcards.sampleID_WiONT]
    output:
        assembly_fa = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly/assembly.fasta",
        assembly_info_txt = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly/assembly_info.txt"
    conda:
       "envs/PacBio_Software_py27_Conda.yml"
    threads: 10
    params:
        Flye_OutputDir_PATH = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly/"
    shell:
        "flye --nano-raw {input.ont_reads_fq} --out-dir {params.Flye_OutputDir_PATH} "
        "--genome-size 4.4m --threads {threads} --asm-coverage 200 --iterations 3"




rule CP_FlyeAssembly_To_I3_Dir: 
    input:
        output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly/assembly.fasta",
    output:
        output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiONT}.flyeassembly.I3.fasta",
    shell: 
        "cp {input} {output}"


### Fitler long read assembly (Flye 3X polished) for only contigs greater than 100kb

rule filterByLength_100kbContigs_FlyeAssembly_I3:
    input:
        PacBio_Flye_Assembly_fa = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiONT}.flyeassembly.I3.fasta",
    output:
        PacBio_Flye_Assembly_Renamed_100KbContigs_FA = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiONT}.flyeassembly.I3.Renamed.100Kb.fasta",
        PacBio_Flye_Assembly_Renamed_100KbContigs_FAI = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiONT}.flyeassembly.I3.Renamed.100Kb.fasta.fai"
    threads: 1
    conda:
        "envs/bioinfo_util_env_V1.yml"
    shell:
        " bioawk -c fastx '{{ print \">{wildcards.sampleID_WiONT}_\"$name \"\\n\" $seq }}' {input.PacBio_Flye_Assembly_fa} "
        " | "
        " bioawk -c fastx '{{ if(length($seq) > 100000) {{ print \">\"$name; print $seq }}}}' > {output.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} \n"
        " samtools faidx {output.PacBio_Flye_Assembly_Renamed_100KbContigs_FA} "



##### Medaka polishing of Flye assembly using ONT reads #####

rule Medaka_Polish_FlyeAsm:
    input:
        ont_reads_fq = lambda wildcards: SampleID_To_ONT_FQ_Dict[wildcards.sampleID_WiONT],
        Flye_Assembly_Renamed_100KbContigs_FA = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Assembly_RenamedAndLengthFiltered/{sampleID_WiONT}.flyeassembly.I3.Renamed.100Kb.fasta",
    output:
        Medaka_Consensus_FA = temp(output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/consensus.fasta"),
        Flye_I3_Medaka_Polished_Renamed_FA = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fasta",
        #ReadToDraft_BAM = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/calls_to_draft.bam",
        #ReadToDraft_BAI = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/calls_to_draft.bam.bai",
    threads: 10
    params:
        Medaka_OutputDir = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka",
        Medaka_Model = "r941_min_high_g303"
    shell:
        "rm {params.Medaka_OutputDir}/calls_to_draft.bam {params.Medaka_OutputDir}/calls_to_draft.bam.bai \n"
        "medaka_consensus -i {input.ont_reads_fq} -d {input.Flye_Assembly_Renamed_100KbContigs_FA} -o {params.Medaka_OutputDir} -t {threads} -m {params.Medaka_Model} \n"
        "cp {output.Medaka_Consensus_FA} {output.Flye_I3_Medaka_Polished_Renamed_FA} "




###################################################################################
######### CIRCLATOR for setting start at DnaA (Assuming Circular genome) ##########
###################################################################################


rule circlator_FixStart_DnaA_UseMedakaPolishedAsm:
    input:
        DnaA_Seq_fa = H37rv_DnaA_FA_PATH,
        Flye_I3_Medaka_Polished_Renamed_FA = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fasta",
    output:
        flye_assembly_FixStart_assembly = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta",
    conda:
        "envs/circulator_151py37_Conda.yml"
    threads: 1
    params:
        circlator_out_prefix = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart"
    shell:
        "circlator fixstart --genes_fa {input.DnaA_Seq_fa} {input.Flye_I3_Medaka_Polished_Renamed_FA} {params.circlator_out_prefix}"

###################################################################################


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

rule bwa_idx_ref_100kbContigs_FlyeAsm_I3M:
    input:
        output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta",
    output:
        output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta.bwt"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "bwa index {input}"


rule bwa_map_IllPE_AlignTo_I3_Assembly:
    input:
        FlyeAsm_I3M_FA = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta",
        FlyeAsm_I3M_FA_BWT_IDX = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta.bwt",
        fq1_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_2_trimmed.fastq",
    output:
        output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.sam"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    params:
        rg=r"@RG\tID:{sampleID_WiONT}\tSM:{sampleID_WiONT}"
    threads: 8
    shell:
        "bwa mem -M -R '{params.rg}' -t {threads} {input.FlyeAsm_I3M_FA} "
        "{input.fq1_trimmed} {input.fq2_trimmed} > {output}"


rule samtools_ViewAndSort_IllPE_AlignTo_I3_Assembly:
    input:
        output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.sam"
    output:
        bam = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.bam",
        bai = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.bam.bai",
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "samtools view -bS {input} | samtools sort - > {output.bam} \n"
        "samtools index {output.bam}"


#####################################
#### PICARD (remove duplicates) #####
#####################################

rule picard_RemoveDup_IllPE_AlnTo_I3M_Asm:
    input:
        IllPE_BwaMEM_BAM = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.bam",
    output:
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.duprem.bam",
        IllPE_BwaMEM_Duprem_BAI = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.duprem.bam.bai",
        IllPE_BwaMEM_Duprem_METRICS = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.duprem.bam.metrics",
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "picard -Xmx5g MarkDuplicates I={input.IllPE_BwaMEM_BAM} O={output.IllPE_BwaMEM_Duprem_BAM} "
        "REMOVE_DUPLICATES=true M={output.IllPE_BwaMEM_Duprem_METRICS} ASSUME_SORT_ORDER=coordinate"
        " \n"
        "samtools index {output.IllPE_BwaMEM_Duprem_BAM}"


rule pilon_IllPE_Polishing_I3M_Asm:
    input:
        FlyeAsm_I3M_FA = output_Dir + "/{sampleID_WiONT}/ONT/Flye_Asm_Medaka/{sampleID_WiONT}.Flye.I3.Medaka.fixstart.fasta",
        IllPE_BwaMEM_Duprem_BAM = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.duprem.bam",
        IllPE_BwaMEM_Duprem_BAI = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/Ill_AlnTo_FlyeAsm_I3M/{sampleID_WiONT}.IllPE.AlnTo.I3MAsm.duprem.bam.bai",
    output:
        pilon_VCF = output_Dir+ "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.vcf",
        I3M_Asm_PilonPolished_FA = output_Dir+ "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta",
        I3M_Asm_PilonPolished_ChangesFile = output_Dir+ "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.changes"
    params:
        Pilon_OutputDir_PATH = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "pilon -Xmx10g --fix snps,indels --genome {input.FlyeAsm_I3M_FA} --bam {input.IllPE_BwaMEM_Duprem_BAM} --output {wildcards.sampleID_WiONT}.Flye.I3MAsm.PilonPolished"
        " --outdir {params.Pilon_OutputDir_PATH} --variant --changes --tracks \n"
        #"samtools faidx {output.pilon_I3_Polished_FA}"



rule samtools_faidx_Flye_I3M_PP_Asm:
    input:
        output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta",
    output:
        output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta.fai",
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    threads: 1
    shell:
        "samtools faidx {input}"



##############################################################################################################
######### CIRCLATOR for setting start at DnaA of Pilon Polished genome (Assuming Circular genome) ############
##############################################################################################################


rule circlator_FixStart_DnaA_Use_I3MPP_Asm:
    input:
        DnaA_Seq_fa = H37rv_DnaA_FA_PATH,
        I3M_Asm_PilonPolished_FA = output_Dir+ "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta",
    output:
        flye_assembly_FixStart_assembly = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fixstart.fasta",
    conda:
        "envs/circulator_151py37_Conda.yml"
    threads: 1
    params:
        circlator_out_prefix = output_Dir+ "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fixstart"
    shell:
        "circlator fixstart --genes_fa {input.DnaA_Seq_fa} {input.I3M_Asm_PilonPolished_FA} {params.circlator_out_prefix}"

##############################################################################################################




###################################################################################
########## Let's run PolyPolish on the Flye-I3-Medaka-PilonPolished Asm ###########
###################################################################################

# Step 1: Create a BWA index for the FlyeI3MPP Asm 

rule bwa_idx_ref_FlyeAsm_I3MPP:
    input:
        output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta",
    output:
        output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta.bwt"
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    shell:
        "bwa index {input}"


rule bwa_mem_AlnAll_ForPolyPolish:
    input:
        FlyeAsm_I3MPP_FA = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta",
        FlyeAsm_I3MPP_FA_BWT_IDX = output_Dir+ "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta.bwt",
        fq1_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_1_trimmed.fastq",
        fq2_trimmed = output_Dir + "/{sampleID_WiONT}/IlluminaWGS/FASTQs_Trimmomatic_Trimming_V2/{sampleID_WiONT}_2_trimmed.fastq",
    output:
        SAM_1 = temp(output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.1.sam"),
        SAM_2 = temp(output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.2.sam"),
    conda:
        "envs/IlluminaPE_Processing_Conda.yml"
    threads: 16
    shell:
        "bwa mem -t 16 -a {input.FlyeAsm_I3MPP_FA} {input.fq1_trimmed} > {output.SAM_1} \n"
        "bwa mem -t 16 -a {input.FlyeAsm_I3MPP_FA} {input.fq2_trimmed} > {output.SAM_2} "


rule polypolish_InsertFilter:
    input:
        SAM_1 = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.1.sam",
        SAM_2 = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.2.sam",
    output:
        SAM_Filtered_1 = temp(output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.Filtered.1.sam"),
        SAM_Filtered_2 = temp(output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.Filtered.2.sam"),
    conda:
        "envs/polypolish_v05.yml"
    shell:
        "polypolish_insert_filter.py --in1 {input.SAM_1} --in2 {input.SAM_2} --out1 {output.SAM_Filtered_1} --out2 {output.SAM_Filtered_2} "



rule polypolish_Polish_FlyeI3MPP_Asm:
    input:
        FlyeAsm_I3MPP_FA = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3M_PilonPolishing/pilon_IllPE_Polishing_I3M_Asm_ChangeSNPsINDELsOnly/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.fasta",
        SAM_Filtered_1 = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.Filtered.1.sam",
        SAM_Filtered_2 = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.bwa.AllAln.Filtered.2.sam",
    output:
        I3MPP_PolyPolished_FA = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.fasta",
        I3MPP_PolyPolished_Debug_TSV = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.Debug.tsv",
        #I3MPP_PolyPolished_Changes_TSV = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.Changes.tsv",
    conda:
        "envs/polypolish_v05.yml"
    shell:
        "polypolish --debug {output.I3MPP_PolyPolished_Debug_TSV} {input.FlyeAsm_I3MPP_FA} {input.SAM_Filtered_1} {input.SAM_Filtered_2} > {output.I3MPP_PolyPolished_FA} \n"





rule circlator_FixStart_DnaA_Use_I3MPPPolyPolish_Asm:
    input:
        DnaA_Seq_fa = H37rv_DnaA_FA_PATH,
        I3MPP_PolyPolished_FA = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.fasta",
    output:
        flye_assembly_FixStart_assembly = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.fixstart.fasta",
    conda:
        "envs/circulator_151py37_Conda.yml"
    threads: 1
    params:
        circlator_out_prefix = output_Dir + "/{sampleID_WiONT}/FlyeAssembly_I3MPP_PolyPolish/{sampleID_WiONT}.Flye.I3MAsm.PilonPolished.PolyPolished.fixstart"
    shell:
        "circlator fixstart --genes_fa {input.DnaA_Seq_fa} {input.I3MPP_PolyPolished_FA} {params.circlator_out_prefix}"





