# 9.Mtb.SRAsms.PGAnalysis.smk
### Snakemake pipeline for performing pan-genome analysis of all Mtb SHORT-READ assemblies (N = 151)



# Read in data regarding input 
input_DataInfo_DF = pd.read_csv( config["inputSampleData_TSV"], sep='\t')


# Create a python list of Sample IDs
input_All_SampleIDs = list( input_DataInfo_DF["SampleID"].values )

SampleID_to_Asm_FA_Dict = dict(input_DataInfo_DF[['SampleID', 'ShortRead_Genome_ASM_PATH']].values)

print("List of input sampleIDs:", len(input_All_SampleIDs), input_All_SampleIDs)



# Define directory with PGAP Annotated Assemblies
PGAP_GenomeAnnoDir = "/n/data1/hms/dbmi/farhat/mm774/Projects/230121_PGAP_AnnoByTBPortals_V1/genomes"

PGAP_SRAsm_GenomeAnnoDir = "/n/data1/hms/dbmi/farhat/mm774/Projects/230313_PGAP_AnnoByTBPortals_SRAsms_V1"




rule all:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.fna", sampleID = input_All_SampleIDs),
        # expand(output_Dir + "/FastANI/FastANI_SR_AsmToH37Rv/{sampleID}.SRAsmToH37Rv.FastANI.txt", sampleID = input_All_SampleIDs),

        expand(output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_SRAsm/{sampleID}/short_summary.specific.corynebacteriales_odb10.{sampleID}.txt", sampleID = input_All_SampleIDs), #["MFS-56"]),

        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates/pangenome.ContentSummary.txt",
        output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_I90_AllIsolates/pangenome.ContentSummary.txt",
        output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_I95_AllIsolates/pangenome.ContentSummary.txt",

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



        # PG Analysis based on PGAP Annotations
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",

        output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.ContentSummary.txt",
        output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_I90_AllIsolates_WiPGAPAnno_V1/pangenome.ContentSummary.txt",
        output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_I95_AllIsolates_WiPGAPAnno_V1/pangenome.ContentSummary.txt",
        






##################################################
############ Genome Assembly QC ##################
##################################################

# Run BUSCO for general QC and completeness evaluation


rule busco_SRAsm:
    input:
        i_Assembly_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        Rv_ShortSum_TXT = output_Dir + "/AsmAnalysis/H37Rv/Assembly_QC/Busco_Mtb_LRAsm/H37Rv/short_summary.specific.corynebacteriales_odb10.H37Rv.txt",
    output:
        ShortSum_TXT = output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_SRAsm/{sampleID}/short_summary.specific.corynebacteriales_odb10.{sampleID}.txt",
    conda:
        "CondaEnvs/Busco_v5_4_4.nobuilds.yml"
    params:
        mode = "genome",
        name = '{sampleID}',
        out_dir = output_Dir + "/AsmAnalysis/{sampleID}/Assembly_QC/Busco_Mtb_SRAsm/",
        tmp_dir = output_Dir + "/Busco_Download_Tmp/",
        lineage="corynebacteriales_odb10",
    shell:
        #"mkdir {params.tmp_dir} \n"
        "mkdir -p {params.out_dir} \n"
        "cd {params.out_dir} \n"
        "busco -i {input.i_Assembly_FA} -o {wildcards.sampleID} --out_path {params.out_dir} --mode {params.mode} "
        " --lineage_dataset {params.lineage} --download_path {params.tmp_dir} --force --offline "




#################################################
############ Genome Annotation ##################
#################################################


rule Bakta_Anno_WiH37Rv_SRAsm:
    input:
        i_Assembly_KrakFilt_FA = lambda wildcards: SampleID_to_Asm_FA_Dict[wildcards.sampleID],
        i_H37Rv_GBK = H37rv_GBK_PATH,
    output:
        Asm_Bakta_Anno_FAA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.faa",
        Asm_Bakta_Anno_FNA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.fna",
        Asm_Bakta_Anno_GBFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gbff",
        Asm_Bakta_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3",
    conda:
        "CondaEnvs/Bakta_1_6_1.nobuilds.yml" 
    threads: 8
    params:
        Bakta_OutputDir_PATH = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/",
        Bakta_DB_Dir = "/n/data1/hms/dbmi/farhat/mm774/References/Bakta_DB_v4/db",
        Bakta_OutPrefix = "{sampleID}.Bakta"
    shell:
        "bakta --db {params.Bakta_DB_Dir} --verbose --output {params.Bakta_OutputDir_PATH} "
        " --prefix {params.Bakta_OutPrefix} --locus-tag {wildcards.sampleID}  "
        " --proteins {input.i_H37Rv_GBK}"
        " --threads {threads} {input.i_Assembly_KrakFilt_FA} "




##################################################
############ Pangenome Analysis ##################
##################################################

rule Panaroo_StdParams_Strict_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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





rule Roary_DefaultParams_AllIsolates_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gff3", sampleID = input_All_SampleIDs),           
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
 































rule create_InputPATHInfo_Bakta_GBFF_SRAsm:
    input:
        Asm_Bakta_Anno_GBFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gbff",
    output:
        Bakta_GBFF_InputPATH_Info_TSV = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gbff.InputInfo.tsv",
    shell:
        'echo "{wildcards.sampleID}\t{input.Asm_Bakta_Anno_GBFF}" > {output.Bakta_GBFF_InputPATH_Info_TSV}'


rule mergeInputGBFF_PATHs_To_TSV_ForPpanggolin_AllIsolates_SRAsm:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_Bakta/{sampleID}.Bakta.gbff.InputInfo.tsv", sampleID = input_All_SampleIDs),       
    output:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Preprocessing/AllGenomes.GBFF.InputPATHs.tsv",
    shell:
        "cat {input} > {output.All_Input_GBFF_PATHs_TSV}"


rule run_Ppanggolin_DefaultParam_AllIsolates_SRAsm:
    input:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Preprocessing/AllGenomes.GBFF.InputPATHs.tsv",
    output:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates/pangenome.h5",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 8
    params:
        Ppanggolin_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates/",
    shell:
        "ppanggolin workflow -f --cpu {threads} --output {params.Ppanggolin_OutputDir} --anno {input.All_Input_GBFF_PATHs_TSV} "

rule get_PangenomeContent_Ppanggolin_AllIsolates_SRAsm:
    input:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates/pangenome.h5",
    output:
        Pangenome_ContentSummary_TXT = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates/pangenome.ContentSummary.txt",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 1
    shell:
        "ppanggolin info -p {input.Pangenome_H5} --content > {output.Pangenome_ContentSummary_TXT}"





##################################################





##################################################################################
##################### Run PGAP Based Pangenome Analysis  #########################
##################################################################################

# /n/data1/hms/dbmi/farhat/mm774/Projects/230313_PGAP_AnnoByTBPortals_SRAsms_V1/R37765/output


rule CP_SR_PGAP_AnnoFiles:
    input:
        i_Anno_GFF = PGAP_SRAsm_GenomeAnnoDir + "/{sampleID}/output/annot.gff", 
        i_Anno_Genome_FA = PGAP_SRAsm_GenomeAnnoDir + "/{sampleID}/output/annot.fna",  
        #i_Anno_FNA = PGAP_SRAsm_GenomeAnnoDir + "/{sampleID}/output/annot.fna", 
        #i_Anno_WiGenomicFASTA_GFF = PGAP_SRAsm_GenomeAnnoDir + "/{sampleID}/output/annot_with_genomic_fasta.gff",        
        #i_Anno_GBK = PGAP_SRAsm_GenomeAnnoDir + "/{sampleID}/output/annot.gbk",  
        #i_Anno_Genome_FA = PGAP_SRAsm_GenomeAnnoDir + "/{sampleID}/output/annot.fna",  
    output:
        Renamed_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.gff", 
        Renamed_Genome_FA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.Genome.fasta",  
    shell: # Rename the GFF & FASTA file to have the SampleID as the contig name
        "cp {input.i_Anno_GFF} {output.Renamed_Anno_GFF} \n"
        ""
        " bioawk -c fastx '{{ print \">\" substr($name, 5) \"\\n\" $seq }}' {input.i_Anno_Genome_FA} > {output.Renamed_Genome_FA}"

        #"cp {input.i_Anno_Genome_FA} {output.Renamed_Genome_FA} \n"
        
        # "awk -F'\\t' '/^#/ {{print; next}} {{OFS=\"\\t\"; $1=\"{wildcards.sampleID}\"; print}}' {input.i_Anno_GFF} >  {output.Renamed_Anno_GFF} \n"
        # " bioawk -c fastx '{{ print \">{wildcards.sampleID}\" \"\\n\" $seq }}' {input.i_Anno_Genome_FA} > {output.Renamed_Genome_FA}"

        #"cp {input.i_Anno_GFF} {output.Renamed_Anno_GFF} \n"
        #""
        #"cp {input.i_Anno_Genome_FA} {output.Renamed_Genome_FA} \n"
        
        # "cp {input.i_Anno_FNA} {output.Renamed_Anno_FNA} \n"
        # "cp {input.i_Anno_FAA} {output.Renamed_Anno_FAA} \n"
        # "cp {input.i_Anno_WiGenomicFASTA_GFF} {output.Renamed_Anno_WiGenomicFASTA_GFF} \n"
        # "cp {input.i_Anno_GBK} {output.Renamed_Anno_GBK} \n"



rule reformat_PGAP_GFF_SRAsm:
    input:
        Renamed_Anno_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.gff", 
        Renamed_Genome_FA = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.Genome.fasta",  
    output:
        Updated_PGAP_GFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", 
    shell:
        "convert_refseq_to_prokka_gff.py -g {input.Renamed_Anno_GFF} -f {input.Renamed_Genome_FA} -o {output.Updated_PGAP_GFF} "



# Run Pangenome analysis tools w/ PGAP Anno Asms

#  PGAP - Panaroo
rule Panaroo_StdParams_Strict_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/summary_statistics.txt",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 4
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes"



rule Panaroo_StdParams_Moderate_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/summary_statistics.txt",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes"



rule Panaroo_StdParams_Sensitive_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/summary_statistics.txt",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes"



# PGAP - Panaroo - Wi Merge Paralogs

##### Run Panaroo but with the '--merge_paralogs' option. 

rule Panaroo_StdParams_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/summary_statistics.txt",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode strict -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/summary_statistics.txt",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 2
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode moderate -t {threads} --remove-invalid-genes --merge_paralogs"


rule Panaroo_StdParams_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       SummaryStats_TXT = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/summary_statistics.txt",
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/Panaroo_1_3_4.nobuilds.yml"
    threads: 8
    params:
        Panaroo_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "mkdir -p {params.Panaroo_OutputDir} \n"
        "panaroo -i {input} -o {params.Panaroo_OutputDir} --clean-mode sensitive -t {threads} --remove-invalid-genes --merge_paralogs"






rule Roary_DefaultParams_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates_WiPGAPAnno_V1/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary -v -r -p {threads} -f {params.Roary_OutputDir} {input} " 
        # -e -n


rule Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
       gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
       #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates_WiPGAPAnno_V1/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary  -v -r -s -p {threads} -f {params.Roary_OutputDir} {input} " # -s means 'dont split paralogs'
        # -e -n
        # -e --mafft



rule Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        #pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates_WiPGAPAnno_V1/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary  -v -r -i 90 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e -n
        # -s means 'dont split paralogs'
        # -i 90 sets the minimumblastp identity to 90%
        # -e --mafft



rule Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_SRAsm:
    input:
       expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff", sampleID = input_All_SampleIDs),    
    output:
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/gene_presence_absence.csv",
        pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/pan_genome_reference.fa",
    conda:
        "CondaEnvs/roary_3_13_WiR.yml"
    threads: 8
    params:
        Roary_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates_WiPGAPAnno_V1/",
    shell:
        'rm -r {params.Roary_OutputDir} \n'
        "roary  -v -r -i 80 -s -p {threads} -f {params.Roary_OutputDir} {input} "
        # -e -n
        # -s means 'dont split paralogs'
        # -i 80 sets the minimumblastp identity to 80%
        # -e --mafft
 






rule create_InputPATHInfo_PGAPAnno_GFF:
    input:
        Asm_PGAP_Anno_GBFF = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.WiDNA.gff",
    output:
        PGAP_GFF_InputPATH_Info_TSV = output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.gff.InputInfo.tsv",
    shell:
        'echo "{wildcards.sampleID}\t{input.Asm_PGAP_Anno_GBFF}" > {output.PGAP_GFF_InputPATH_Info_TSV}'


rule mergeInputGBFF_PATHs_To_TSV_ForPpanggolin_AllIsolates_WiPGAPAnno:
    input:
        expand(output_Dir + "/AsmAnalysis/{sampleID}/GenomeAnnotation/{sampleID}_SR_Asm_PGAP_V1/{sampleID}.PGAP.gff.InputInfo.tsv", sampleID = input_All_SampleIDs),       
    output:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Preprocessing_PGAP/AllGenomes.GFF.InputPATHs.tsv",
    shell:
        "cat {input} > {output.All_Input_GBFF_PATHs_TSV}"


rule run_Ppanggolin_DefaultParam_AllIsolates_WiPGAPAnno_SRAsm:
    input:
        All_Input_GBFF_PATHs_TSV = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Preprocessing_PGAP/AllGenomes.GFF.InputPATHs.tsv",
    output:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.h5",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 8
    params:
        Ppanggolin_OutputDir = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/",
    shell:
        "ppanggolin workflow -f --cpu {threads} --output {params.Ppanggolin_OutputDir} --anno {input.All_Input_GBFF_PATHs_TSV} "

rule get_PangenomeContent_Ppanggolin_AllIsolates_WiPGAPAnno_SRAsm:
    input:
        Pangenome_H5 = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.h5",
    output:
        Pangenome_ContentSummary_TXT = output_Dir + "/PanGenome_Analysis/SR_Ppanggolin_Default_AllIsolates_WiPGAPAnno_V1/pangenome.ContentSummary.txt",
    conda:
        "CondaEnvs/ppanggolin_1_2_74.yml"
    threads: 1
    shell:
        "ppanggolin info -p {input.Pangenome_H5} --content > {output.Pangenome_ContentSummary_TXT}"











###### Run panqc analysis steps #####



rule Panaroo_Strict_MergeParalogs_AvA_KmerComparison_SRAsm:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "

rule Panaroo_Moderate_MergeParalogs_AvA_KmerComparison_SRAsm:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp,
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "

rule Panaroo_Sensitive_MergeParalogs_AvA_KmerComparison_SRAsm:
    input:
       pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "




rule Roary_DefaultParams_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "

rule Roary_NoSplitParalogs_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "


rule Roary_NoSplitParalogs_I90_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "



rule Roary_NoSplitParalogs_I80_AvA_KmerComparison:
    input:
       pangenome_ref_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PG_AvA_MaxJC_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.KmerComparison.AllVsAll.MaxJC.tsv",   
    params:
       k = 31 # k-mer size in bp
    shell:
        "time panqc ava -k {params.k} -i {input} -o {output} "





###### Run GeneSeqCheck step for all Panaroo & Roary analyses ######

rule create_SampleToAsmFA_TSV_SRAsm:
    output:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
    run:
        # Output as tsv
        input_DataInfo_SRAsm_DF = input_DataInfo_DF[['SampleID', 'ShortRead_Genome_ASM_PATH']].copy()
        input_DataInfo_SRAsm_DF.columns = ['SampleID', 'Genome_ASM_PATH']
        input_DataInfo_SRAsm_DF.to_csv(output.asmfa_tsv, sep="\t", index = False)


rule Panaroo_Strict_MergeParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Strict_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Panaroo_Moderate_MergeParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Moderate_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Panaroo_Sensitive_MergeParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Panaroo_Sensitive_MergeParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Roary_DefaultParams_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_Default_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",  
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Roary_NoSplitParalogs_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


rule Roary_NoSplitParalogs_I90_panqc_AsmGeneSeqCheck:
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I90_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "



rule Roary_NoSplitParalogs_I80_panqc_AsmGeneSeqCheck: 
    input:
        asmfa_tsv = output_Dir + "/PanGenome_Analysis/InputAssembly.ShortRead.FA.PATHs.tsv",
        gene_presence_absence_CSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.csv",
        pan_genome_reference_FA = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/pan_genome_reference.fa",
    output:
       PresAbs_AsmSeqCheck_TSV = output_Dir + "/PanGenome_Analysis/SR_Roary_NoSplitParalogs_I80_AllIsolates/gene_presence_absence.AsmGeneSeqChk.tsv",
    shell:
        "time panqc asmseqcheck -a {input.asmfa_tsv} "
        " -m {input.gene_presence_absence_CSV} "
        " -r {input.pan_genome_reference_FA} "
        " -o {output} "


#####################################################################









