# Minigraph analysis of H37Rv (Reference) & 151 Mtb genomes (complete hybrid assemblies)


### Import Statements ###
import pandas as pd


### Define PATHs to files defined in thoe config file ###
refGenome_FA_PATH = config["RefGenome_FA_PATH"]

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
        output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.V1.gfa",
        expand(output_Dir + "/AsmAnalysis/{sampleID}/Minigraph_CallSVs/{sampleID}.Minigraph.CallSVs.ToH37RvAnd158CI.bed", sampleID = input_All_SampleIDs),
        output_Dir + "/AsmAnalysis/H37Rv/Minigraph_CallSVs/H37Rv.Minigraph.CallSVs.ToH37RvAnd158CI.bed",
        
        output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.MergedSV.Info.PASTE.tsv",

        output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.MergedSV.Info.svvcf",
        output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.MergedSV.Info.tsv",



############################################################################
############### Minigraph: Pangenome construction & Analysis ###############
############################################################################

rule create_TMP_Input_Asm_FA:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
    output:
        tmp_Asm_Renamed_FA = temp(output_Dir + "/Minigraph/tmp/{sampleID}.LR.Asm.Renamed.fasta")
    conda:
        "CondaEnvs/Minigraph_v0_19.yml"
    shell:
        " bioawk -c fastx '{{ print \">{wildcards.sampleID}\" \"\\n\" $seq }}' {input} > {output}"


# sampleID = ["DNA086", "M0011368_9"],#

rule run_Minigraph_WiH37RvRef:
    input:
        refGenome_FA_PATH,
        expand( output_Dir + "/Minigraph/tmp/{sampleID}.LR.Asm.Renamed.fasta", sampleID = input_All_SampleIDs),
    output:
        MG_PG_GFA = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.V1.gfa",
        MG_PG_Bubble_SV_BED = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.V1.Bubble.SV.bed",
        MG_PG_Stable_FA = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.V1.Stable.fa",
    #conda:
    #    "CondaEnvs/Minigraph_v019_gfatools_05.yml"
    threads: 4
    shell:
        "time minigraph -cxggs -t {threads} {input} > {output.MG_PG_GFA} \n"
        "gfatools bubble {output.MG_PG_GFA} > {output.MG_PG_Bubble_SV_BED} \n"
        "gfatools gfa2fa -s {output.MG_PG_GFA}  > {output.MG_PG_Stable_FA} \n"

rule minigraph_CallSVs:
    input:
        i_Assembly_FA = output_Dir + "/Minigraph/tmp/{sampleID}.LR.Asm.Renamed.fasta",
        MG_PG_GFA = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.V1.gfa"
    output:
        CallSV_BED = output_Dir + "/AsmAnalysis/{sampleID}/Minigraph_CallSVs/{sampleID}.Minigraph.CallSVs.ToH37RvAnd158CI.bed",
    #conda:
    #    "CondaEnvs/Minigraph_v019_gfatools_05.yml"
    threads: 1
    shell: 
        "time minigraph -cxasm --call {input.MG_PG_GFA} {input.i_Assembly_FA} > {output.CallSV_BED} "


rule minigraph_CallSVs_H37Rv:
    input:
        i_Assembly_FA = refGenome_FA_PATH,
        MG_PG_GFA = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.V1.gfa"
    output:
        CallSV_BED = output_Dir + "/AsmAnalysis/H37Rv/Minigraph_CallSVs/H37Rv.Minigraph.CallSVs.ToH37RvAnd158CI.bed",
    #conda: "CondaEnvs/Minigraph_v019_gfatools_05.yml"
    threads: 1
    shell: 
        "time minigraph -cxasm --call {input.MG_PG_GFA} {input.i_Assembly_FA} > {output.CallSV_BED} "




rule paste_Minigraph_SV_Calls:
    input:
        output_Dir + "/AsmAnalysis/H37Rv/Minigraph_CallSVs/H37Rv.Minigraph.CallSVs.ToH37RvAnd158CI.bed",
        expand(output_Dir + "/AsmAnalysis/{sampleID}/Minigraph_CallSVs/{sampleID}.Minigraph.CallSVs.ToH37RvAnd158CI.bed", sampleID = input_All_SampleIDs),
    output:
        paste_SV_BEDs_TSV = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.MergedSV.Info.PASTE.tsv",
    shell: 
        "paste {input} > {output}"

rule merge_Minigraph_SV_Calls:
    input:
        output_Dir + "/AsmAnalysis/H37Rv/Minigraph_CallSVs/H37Rv.Minigraph.CallSVs.ToH37RvAnd158CI.bed",
        expand(output_Dir + "/AsmAnalysis/{sampleID}/Minigraph_CallSVs/{sampleID}.Minigraph.CallSVs.ToH37RvAnd158CI.bed", sampleID = input_All_SampleIDs),
    output:
        MG_PG_MergedSV_Info_TSV = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.MergedSV.Info.tsv",
        MG_PG_MergedSV_Info_SVVCF = output_Dir + "/Minigraph/Minigraph_H37rv_Vs_158CI.MergedSV.Info.svvcf",
    shell: 
        "paste {input} | mgutils.js merge - > {output.MG_PG_MergedSV_Info_TSV} \n"
        "mgutils.js merge2vcf < {output.MG_PG_MergedSV_Info_TSV} > {output.MG_PG_MergedSV_Info_SVVCF} "











