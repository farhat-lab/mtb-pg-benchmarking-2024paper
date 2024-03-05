
# PMP_SM - Preprocessing for merging PacBio FQs (Any # of FQs to single merged FQ for one isolate)


# This is a test run of the merging approach for FQs


### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]
refGenome_GFF_PATH = config["RefGenome_GFF_PATH"]

H37rv_DnaA_FA_PATH = config["H37rv_DnaA_FA_PATH"]
H37rv_GBK_PATH = config["H37rv_GBK_PATH"]

Mcanetti_VCF_PATH = config["Mcanettii_VCF"]
Mcanetti_SNP_Positions_TSV_PATH = config["Mcanettii_SNP_Positions_TSV"]


# Define PATH of main OUTPUT directory
output_Dir = config["output_dir"]



# Define analysis name to use a prefix for output files
#AnalysisName = config["analysis_name"]

# Define window sizes for nucleotide diversity calculation
NucDiv_WindowSizes_bp = ['1000']



# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')


# Create a python list of Sample IDs
input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )

SampleID_to_Asm_FA_Dict = dict(input_DataInfo_DF[['SampleID', 'Genome_ASM_PATH']].values)

print("List of input sampleIDs:", len(input_All_SampleIDs), input_All_SampleIDs)


rule all:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.var.tsv", sampleID = input_All_SampleIDs),

        output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv",
        output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",

        output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.min-100.fasta",
        output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.fasta",
        output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.min-100.fasta",


        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_NoFilt/MM2.mpileup.call.Merged.SNVs.min-100.fasttree.newick",   
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasttree.newick",
        
        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq.treefile",
        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.iq.treefile",
        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.iq.treefile",
        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_RLCandLowPmapMask/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.iq.treefile",

        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/IQtree.10AmbThresh.PLCMask.MidRoot.WiNodeNames.newick",
        output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/IQtree.10AmbThresh.NoMask.MidRoot.WiNodeNames.newick",



#########################################################
#########################################################
## Analysis Versus H37rv (call variants against H37Rv) ##
#########################################################
#########################################################



##############################################################################
############ Minimap2: Assembly To H37rv Alignment & Variant Calling #########
##############################################################################

rule MM2_AsmToH37rv:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        H37rv_FA = refGenome_FA_PATH,
    output:
        MM2_AsmToH37rv_SAM = temp(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.sam"),
        MM2_AsmToH37rv_BAM = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.bam",
        MM2_AsmToH37rv_bai = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.bam.bai",
        MM2_AsmToH37rv_paftools_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.paftools.vcf",
        MM2_AsmToH37rv_VarTSV = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.var.tsv",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    threads: 1
    params:
        MM2_MinAlnLen_ForCoverage = 1000,
        MM2_MinAlnLen_ForVariantCalling = 1000,
    shell: 
        "minimap2 -ax asm10 --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^@/ || ($5 == 60)' > {output.MM2_AsmToH37rv_SAM} \n"
        "samtools view -bS {output.MM2_AsmToH37rv_SAM} | samtools sort - > {output.MM2_AsmToH37rv_BAM} \n"
        "samtools index {output.MM2_AsmToH37rv_BAM} \n"
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^R/ || ($12 == 60)' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} -f {input.H37rv_FA} - > {output.MM2_AsmToH37rv_paftools_VCF} \n"                         
        "minimap2 -cx asm10 --cs {input.H37rv_FA} {input.i_Assembly_FA} | awk '$1 ~ /^R/ || ($12 == 60)' | sort -k6,6 -k8,8n | paftools.js call -s {wildcards.sampleID} -L {params.MM2_MinAlnLen_ForVariantCalling} -l {params.MM2_MinAlnLen_ForCoverage} - > {output.MM2_AsmToH37rv_VarTSV} \n"




##########################################################################
############ Subset SNPs from Mpileup for Phylogeny building #############

rule bcftools_mpileup_MM2_AsmToH37rv:
    input:
        MM2_AsmToH37rv_BAM = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37Rv.bam",
        H37rv_FA = refGenome_FA_PATH,     
    output:
        MM2_AsmToH37rv_BAM_mpileup_out = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.txt.vcf",
        MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.KeepAllPositions.vcf"
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    threads: 1
    shell: # Run bcftools mpileup to summarize coverage and basepair outputs from the Minimap2 alignment (BAM)
        "bcftools mpileup -f {input.H37rv_FA} {input.MM2_AsmToH37rv_BAM} > {output.MM2_AsmToH37rv_BAM_mpileup_out} \n"
        "bcftools mpileup -f {input.H37rv_FA} {input.MM2_AsmToH37rv_BAM} | bcftools call -c -o {output.MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF}"



##### Merge SNPs for Phylogeny generation (Based on Flye Asm) #####

rule getAll_SNPpositions_mpileup_MM2_AsmToH37rv:
    input:
        MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.KeepAllPositions.vcf",
    output:
        MM2_mpileup_VCF_AllSNPpositions_TSV = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.AllSNPpositions.tsv",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        "bcftools view --types snps {input.MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF} | cut -f 1,2 | grep -v '#' > {output.MM2_mpileup_VCF_AllSNPpositions_TSV} \n"


rule combineAll_SNPpositions_mpileup_MM2_AsmToH37rv:
   input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.AllSNPpositions.tsv", sampleID = input_All_SampleIDs),           
   output:
       AllSample_AllSNPpositions_MM2_mpileup_TSV = output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv",
   shell:
       "cat {input} | sort -k 2n | uniq > {output.AllSample_AllSNPpositions_MM2_mpileup_TSV}"


rule Filter_MM2_mpileup_VarCalling_To_OnlySNPpositionsInUnionOfAllSNPs:
    input:
        MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.bam.mpileup.call.KeepAllPositions.vcf",
        AllSample_AllSNPpositions_MM2_mpileup_TSV = output_Dir + "/Asm_MergeSNPs_mpileup/AllSNPpositions.MM2.mpileup.SNPs.Union.AllSamples.tsv"
    output:
        MM2_AsmToRef_AllPositions_BCF_GZ = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.AllPositions.bcf.gz",
        MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.mpileup.call.SNPs.Union.AllSamples.bcf.gz",
        SampleID_TXT = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.name.txt",
        MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.mm2.AsmToH37rv.mpileup.call.SNPs.Union.AllSamples.Renamed.bcf.gz",
        BCF_Renamed_PATH_TXT = output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.BCF_RenamedPATH.txt",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        "bcftools view {input.MM2_AsmToH37rv_BAM_mpileup_call_KeepAllPositions_VCF} -O b -o {output.MM2_AsmToRef_AllPositions_BCF_GZ} \n"
        
        "bcftools index {output.MM2_AsmToRef_AllPositions_BCF_GZ} \n"

        "bcftools view {output.MM2_AsmToRef_AllPositions_BCF_GZ} "
        " -R {input.AllSample_AllSNPpositions_MM2_mpileup_TSV} -e 'DP!=1' -O b -o {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples} \n"

        "bcftools index {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples} \n"

        "echo {wildcards.sampleID} > {output.SampleID_TXT} \n"
        
        "bcftools reheader -s {output.SampleID_TXT} {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples} -o {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed} \n "
        
        "bcftools index {output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed} \n"
        
        "echo '{output.MM2_AsmToRef_SNPs_BCF_GZ_SNPsInAllSamples_Renamed}' > {output.BCF_Renamed_PATH_TXT} "



rule merge_BCFs_Renamed_PATHs_To_TXT:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/VariantCallingVersusH37Rv/MM2_AsmToH37rv/{sampleID}.BCF_RenamedPATH.txt", sampleID = input_All_SampleIDs), 
    output:
        output_Dir + "/Asm_MergeSNPs_mpileup/ListOfAll_PATHs_BCFs_Renamed.txt"
    shell:
        "cat {input} > {output}"



rule merge_All_BCFs_Renamed_To_VCF:
    input:
        output_Dir + "/Asm_MergeSNPs_mpileup/ListOfAll_PATHs_BCFs_Renamed.txt"
    output:
        Merged_VCF = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.vcf",
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
        Merged_VCF_SNVs_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.MergedSNPs.snp.positions",
        Merged_VCF_SNVs_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf.gz",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        ' bcftools merge -i "-" -l {input} -o {output.Merged_VCF} -O v \n' # Does NOT fill AMB positions, ambigous calls stay ambigous
        ""
        " bcftools view {output.Merged_VCF} --types snps > {output.Merged_VCF_SNVs} \n"

        'grep -v "#" {output.Merged_VCF_SNVs} | cut -f 2  > {output.Merged_VCF_SNVs_POS} \n'
        ""
        " bgzip -c {output.Merged_VCF_SNVs} > {output.Merged_VCF_SNVs_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_GZ} "



rule convert_MergedVCF_To_FASTA_ALN:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
    output:
        MergedSNVs_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.min-100.fasta",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    params:
        min_Supporting = -100 
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"



### Filtering of Merged-VCF ### 


### Filter by AMB threshold (Maximum % of AMB allowed a position)

rule filter_SNVs_10AmbThresh_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf.gz"
    output:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
        Merged_VCF_SNVs_10AmbFilt_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf.gz",
        Merged_VCF_SNVs_10AmbFilt_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.positions",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        'bcftools view {input.Merged_VCF_SNVs_GZ} --types snps -e "F_MISSING > 0.10"  > {output.Merged_VCF_SNVs_10AmbFilt} \n'

        'grep -v "#" {output.Merged_VCF_SNVs_10AmbFilt} | cut -f 2  > {output.Merged_VCF_SNVs_10AmbFilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_10AmbFilt} > {output.Merged_VCF_SNVs_10AmbFilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_10AmbFilt_GZ} "


########################################################################



rule filter_SNVs_RLC_Regions_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
        RLC_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.H37Rv.bed"
    output:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.vcf",
        Merged_VCF_SNVs_RLCfilt_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.vcf.gz",
        Merged_VCF_SNVs_RLCfilt_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.positions",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        "bedtools intersect -header -v -a {input.Merged_VCF_SNVs} -b {input.RLC_regions_BED} -wa > {output.Merged_VCF_SNVs_RLCfilt} \n"
        'grep -v "#" {output.Merged_VCF_SNVs_RLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_RLCfilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_RLCfilt} > {output.Merged_VCF_SNVs_RLCfilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_RLCfilt_GZ} "



rule filter_SNVs_PLC_Regions_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.vcf",
        PLC_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/PLC_Regions.H37Rv.bed"
    output:
        Merged_VCF_SNVs_PLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.vcf",
        Merged_VCF_SNVs_PLCfilt_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.vcf.gz",
        Merged_VCF_SNVs_PLCfilt_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.positions",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        "bedtools intersect -header -v -a {input.Merged_VCF_SNVs} -b {input.PLC_regions_BED} -wa > {output.Merged_VCF_SNVs_PLCfilt} \n"
        'grep -v "#" {output.Merged_VCF_SNVs_PLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_PLCfilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_PLCfilt} > {output.Merged_VCF_SNVs_PLCfilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_PLCfilt_GZ} "


rule filter_SNVs_RLCandLowPMap_Regions_MergeSNVs_mpileup:
    input:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
        RLCandLowPmap_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed"
    output:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.vcf",
        Merged_VCF_SNVs_RLCfilt_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.vcf.gz",
        Merged_VCF_SNVs_RLCfilt_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.positions",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        "bedtools intersect -header -v -a {input.Merged_VCF_SNVs_10AmbFilt} -b {input.RLCandLowPmap_regions_BED} -wa > {output.Merged_VCF_SNVs_RLCfilt} \n"
        'grep -v "#" {output.Merged_VCF_SNVs_RLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_RLCfilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_RLCfilt} > {output.Merged_VCF_SNVs_RLCfilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_RLCfilt_GZ} "



rule filter_SNVs_PLC_Regions_MergeSNVs_mpileup_10AmbFilt:
    input:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
        PLC_regions_BED = "References/Mtb_H37Rv_MaskingSchemes/PLC_Regions.H37Rv.bed"
    output:
        Merged_VCF_SNVs_10AmbFilt_PLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.vcf",
        Merged_VCF_SNVs_10AmbFilt_PLCfilt_GZ = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.vcf.gz",
        Merged_VCF_SNVs_10AmbFilt_PLCfilt_POS = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.positions",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    shell:
        "bedtools intersect -header -v -a {input.Merged_VCF_SNVs_10AmbFilt} -b {input.PLC_regions_BED} -wa > {output.Merged_VCF_SNVs_10AmbFilt_PLCfilt} \n"
        'grep -v "#" {output.Merged_VCF_SNVs_10AmbFilt_PLCfilt} | cut -f 2  > {output.Merged_VCF_SNVs_10AmbFilt_PLCfilt_POS} \n'

        " bgzip -c {output.Merged_VCF_SNVs_10AmbFilt_PLCfilt} > {output.Merged_VCF_SNVs_10AmbFilt_PLCfilt_GZ} \n"
        " tabix {output.Merged_VCF_SNVs_10AmbFilt_PLCfilt_GZ} "




rule convert_MergedVCF_10AmbThresh_PLCFilt_To_FASTA_ALN:
    input:
        Merged_VCF_SNVs_10AmbFilt_PLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.vcf",
    output:
        MergedSNVs_10AmbFilt_PLCfilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.fasta",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    params:
        min_Supporting = -100 
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"




rule convert_MergedVCF_10AmbThresh_To_FASTA_ALN:
    input:
        Merged_VCF_SNVs_10AmbFilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.vcf",
    output:
        MergedSNVs_10AmbFilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasta",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    params:
        min_Supporting = -100 
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"




rule convert_MergedVCF_To_FASTA_ALN_RLC_removed:
    input:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.vcf",
    output:
        MergedSNVs_RLCfilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.RLCMask.min-100.fasta",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    params:
        min_Supporting = -100 
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"



rule convert_MergedVCF_To_FASTA_ALN_PLC_removed:
    input:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.vcf",
    output:
        MergedSNVs_PLCfilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.fasta",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    params:
        min_Supporting = -100 
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"


rule convert_MergedVCF_To_FASTA_ALN_RLCandLowPmap_removed:
    input:
        Merged_VCF_SNVs_RLCfilt = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.vcf",
    output:
        MergedSNVs_RLCandPmap_filt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.fasta",
    conda:
        "CondaEnvs/mm2_v2_4_WiUtilities.yml"
    params:
        min_Supporting = -100 
    shell:
        "Scripts/vcf2phylip/vcf2phylip.py -i {input} -f -m {params.min_Supporting} \n"







######## Phylogeny Building ########

rule fasttree_GTR_from_MergedSNPs:
    input:
        MergedSNVs_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.min-100.fasta",
    output:
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_NoFilt/MM2.mpileup.call.Merged.SNVs.min-100.fasttree.newick"   
    conda:
        "CondaEnvs/Gubbins_v3_2.yml"
    shell:
        "time FastTree -nt -gtr {input} > {output}"


rule fasttree_GTR_from_MergedSNVs_10AmbFilt:
    input:
        MergedSNVs_10AmbFilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasta",
    output:
        output_Dir + "/Phylogenies/fasttree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasttree.newick",
    conda:
        "CondaEnvs/Gubbins_v3_2.yml"
    shell:
        "time FastTree -nt -gtr {input} > {output}"

####################################



rule IQtree_GTR_from_MergedSNVs_10AmbFilt:
    input:
        MergedSNVs_10AmbFilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.fasta",
    output:
        MergedSNVs_10AmbFilt_SNPSites_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq.treefile",
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq",
        bb = 10000
    shell:
        "snp-sites {input.MergedSNVs_10AmbFilt_FA} > {output.MergedSNVs_10AmbFilt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_10AmbFilt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo
        #"time iqtree -m GTR -nt 1 -bb {params.bb} -s {input} -pre {params.iqtree_outprefix} " # -m # -redo


rule Update_IQtree_Newick_10AmbNoMask_Filt:
    input:
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.NoMask.min-100.iq.treefile",
    output:
        TREE_WiNodeNames = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_10AmbFilt/IQtree.10AmbThresh.NoMask.MidRoot.WiNodeNames.newick",
    shell:
        "Scripts/IQTree.UpdateNodeNames.V1.py -i {input} -o {output}"




rule IQtree_GTR_from_10AmbFilt_MergedSNVs_PLCFilt:
    input:
        MergedSNVs_10AmbFilt_PLCfilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.fasta",
    output:
        MergedSNVs_10AmbFilt_PLCfilt_SNPSites_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.iq.treefile",
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.iq",
        bb = 10000
    shell:
        "snp-sites {input.MergedSNVs_10AmbFilt_PLCfilt_FA} > {output.MergedSNVs_10AmbFilt_PLCfilt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_10AmbFilt_PLCfilt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo
        #"time iqtree -m GTR -nt 1 -bb {params.bb} -s {input} -pre {params.iqtree_outprefix} " # -m # -redo

rule Update_IQtree_Newick_10AmbPLCMask_Filt:
    input:
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.10AmbThresh.PLCMask.min-100.iq.treefile",
    output:
        TREE_WiNodeNames = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/IQtree.10AmbThresh.PLCMask.MidRoot.WiNodeNames.newick",
    shell:
        "Scripts/IQTree.UpdateNodeNames.V1.py -i {input} -o {output}"



rule IQtree_GTR_from_MergedSNVs_PLCFilt:
    input:
        MergedSNVs_PLCfilt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.fasta",
    output:
        MergedSNVs_PLCfilt_SNPSites_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.iq.treefile",
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/MM2.mpileup.call.Merged.SNVs.PLCMask.min-100.iq",
        bb = 10000
    shell:
        "snp-sites {input.MergedSNVs_PLCfilt_FA} > {output.MergedSNVs_PLCfilt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_PLCfilt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo
        #"time iqtree -m GTR -nt 1 -bb {params.bb} -s {input} -pre {params.iqtree_outprefix} " # -m # -redo






rule IQtree_GTR_from_MergedSNVs_10AmbFilt_RLCandLowPmap_Filt:
    input:
        MergedSNVs_RLCandPmap_filt_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.fasta",
    output:
        MergedSNVs_RLCandPmap_filt_SNPSites_FA = output_Dir + "/Asm_MergeSNPs_mpileup/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.VarSites.fasta",
        TREE = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_RLCandLowPmapMask/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.iq.treefile",
    params:
        iqtree_outprefix = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_RLCandLowPmapMask/MM2.mpileup.call.Merged.SNVs.10AmbThresh.RLCandLowPmapMask.min-100.iq",
        bb = 10000
    shell:
        "snp-sites {input.MergedSNVs_RLCandPmap_filt_FA} > {output.MergedSNVs_RLCandPmap_filt_SNPSites_FA} \n"
        "time iqtree -m GTR+ASC -nt 1 -bb {params.bb} -s {output.MergedSNVs_RLCandPmap_filt_SNPSites_FA} -pre {params.iqtree_outprefix} " # -m # -redo
        #"time iqtree -m GTR -nt 1 -bb {params.bb} -s {input} -pre {params.iqtree_outprefix} " # -m # -redo





##### Run Panstripe Phylo Analysis w/ GAIN + LOSS measure #####

# Panstripe rule in progress here!!!!

rule RunPanstripe_Panaroo_Strict_MergeParalogs:
    input:
       PA_Rtab = output_Dir + "/PanGenome_Analysis/Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.Rtab",
       in_tree = output_Dir + "/Phylogenies/iqtree_mpileupSNVs_PLCFilt/IQtree.10AmbThresh.PLCMask.MidRoot.WiNodeNames.newick",
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















