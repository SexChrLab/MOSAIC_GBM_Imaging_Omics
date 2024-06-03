import os

configfile: "RNAFastqFilePaths_config.rdgroups.json"

# Tool paths:
fastqc_path = "fastqc"
multiqc_path = "multiqc"
samtools_path = "samtools"
bbduksh_path = "/opt/bbmap/bbduk.sh"
hisat_path = "/opt/hisat2-2.2.1/hisat2"
bamtools_path = "bamtools"
picard_path = "picard"
featureCounts_path = "featureCounts"

adapter_path = "/references/Illumina_adapter_sequence.fa"

T2Tv2_Transcriptome_Index_HISAT_Path_male = "/references/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded_hisat"
T2Tv2_Transcriptome_Index_HISAT_Path_female = "/references/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded_hisat"
T2Tv2_annotation_gff_path = "/references/CHM13.v2.0.gff3.gz"

# Directory to put all the results
output_directory = "/data/hawkinsdaarud/RNA/fastqc/"
#output_directory = "/data/Test_variant_callers/RNAseq/"

rule all:
   input:
        ##expand(output_directory+"fastqc_results_rna/{sample_name}_R1_fastqc.html", sample_name = config["all_rna_samples"]),
        ##expand(output_directory+"fastqc_results_rna/{sample_name}_R2_fastqc.html", sample_name = config["all_rna_samples"]),
        #expand(output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R1.fastq.gz", sample_name = config["all_rna_samples"]),
        #expand(output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R2.fastq.gz", sample_name = config["all_rna_samples"]),
        ##expand(output_directory+"fastqc_results_trimmed_rna/{sample_name}_trimmed_R1_fastqc.html", sample_name = config["all_rna_samples"]),
        ##expand(output_directory+"fastqc_results_trimmed_rna/{sample_name}_trimmed_R2_fastqc.html", sample_name = config["all_rna_samples"])
        #expand(output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_XY.sam", male_sample = config["male_samples"]),
        #expand(output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_XX.sam", female_sample = config["female_samples"]),
        #expand(output_directory+"feature_counts_rna/sccREF/{male_sample}_HISAT_geneCounts_XY.txt", male_sample = config["male_samples"]),
        #expand(output_directory+"feature_counts_rna/sccREF/{female_sample}_HISAT_geneCounts_XX.txt", female_sample = config["female_samples"]),
        expand(output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R1.fastq.gz", sample_name = config["female_samples1"]),
        #expand(output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_XY.sam", male_sample= config["male_samples"]),
        #expand(output_directory+"feature_counts_rna/sccREF/{male_sample}_HISAT_geneCounts_XY.txt", male_sample = config["male_samples"])
        expand(output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_XX.sam", female_sample = config["female_samples1"])
        expand(output_directory+"feature_counts_rna/sccREF/{female_sample}_HISAT_geneCounts_XX.txt", female_sample = config["female_samples1"])

#rule fastqc_analysis_untrimmed:
#    input:
#        original_R1 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq_1"],
#        original_R2 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq_2"]
#    output:
#        F1 = output_directory+"fastqc_results_rna/{sample_name}_R1_fastqc.html",
#        F2 = output_directory+"fastqc_results_rna/{sample_name}_R2_fastqc.html",
#    params: 
#        path = fastqc_path,
#        outputdirectory = output_directory+"fastqc_results_rna/"
#    shell:
#        """
#        {params.path} -o {params.outputdirectory} {input.original_R1}, 
#        {params.path} -o {params.outputdirectory} {input.original_R2} 
#        """

rule trim_adapters_paired_bbduk:
    input:
        original_R1 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq_1"],
        original_R2 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq_2"]
    output:
        out_fq_1 = output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R1.fastq.gz",
        out_fq_2 = output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R2.fastq.gz"
    params:
        bbduksh = bbduksh_path,
        adapter = adapter_path
    shell:
        "{params.bbduksh} -Xmx3g in1={input.original_R1} in2={input.original_R2} "
        "out1={output.out_fq_1} out2={output.out_fq_2} "
        "ref={params.adapter} "
        "trimq=20 minlen=30"

#rule fastqc_analysis_trimmed:
#    input:
#        R1 = output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R1.fastq.gz",
#        R2 = output_directory+"trimmed_fastqs_rna/{sample_name}_trimmed_R2.fastq.gz"
#    output:
#        F1 = output_directory+"fastqc_results_trimmed_rna/{sample_name}_trimmed_R1_fastqc.html",
#        F2 = output_directory+"fastqc_results_trimmed_rna/{sample_name}_trimmed_R2_fastqc.html"
#    params: 
#        path = fastqc_path,
#        outputdirectory = output_directory+"fastqc_results_trimmed_rna/"
#    shell:
#        """
#        {params.path} -o {params.outputdirectory} {input.R1}, 
#        {params.path} -o {params.outputdirectory} {input.R2} 
#        """
#
rule HISAT_paired_males:
    input:
        Trimmed_FASTQ1 = output_directory+"trimmed_fastqs_rna/{male_sample}_trimmed_R1.fastq.gz",
        Trimmed_FASTQ2 = output_directory+"trimmed_fastqs_rna/{male_sample}_trimmed_R2.fastq.gz"
    output:
        out_1 = output_directory+ "processed_bams/rna/"+"{male_sample}_HISAT_pair_trim_XY.sam"
    params:
        HISAT_Index_male = T2Tv2_Transcriptome_Index_HISAT_Path_male,
        hisat = hisat_path
    shell:
        "{params.hisat} -q --rna-strandness RF -p 8 -x {params.HISAT_Index_male} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"

rule HISAT_paired_females:
    input:
        Trimmed_FASTQ1 = output_directory+"trimmed_fastqs_rna/{female_sample}_trimmed_R1.fastq.gz",
        Trimmed_FASTQ2 = output_directory+"trimmed_fastqs_rna/{female_sample}_trimmed_R2.fastq.gz"
    output:
        out_1 = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_XX.sam"
    params:
        HISAT_Index_female = T2Tv2_Transcriptome_Index_HISAT_Path_female,
        hisat = hisat_path
    shell:
        "{params.hisat} -q --rna-strandness RF -p 8 -x {params.HISAT_Index_female} -1 {input.Trimmed_FASTQ1} -2 {input.Trimmed_FASTQ2} -S {output.out_1}"


rule samtools_view_males:
    input:
        SAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_XY.sam"
    output:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_XY.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.SAM} > {output.BAM}"

rule bam_sort_males:    
    input:
        IN_BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_XY.bam"
    output:
        sort_BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_XY.bam"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_males:
    input:
        sort_BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_XY.bam"
    output:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam",
        metrics = output_directory+"stats/{male_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_males:
    input:
        Read_BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_XY.bam"
    output:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    params:
        id = lambda wildcards: config[wildcards.male_sample]["ID"],
        sm = lambda wildcards: config[wildcards.male_sample]["SM"],
        lb = lambda wildcards: config[wildcards.male_sample]["LB"],
        pu = lambda wildcards: config[wildcards.male_sample]["PU"],
        pl = lambda wildcards: config[wildcards.male_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_males:
    input:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        BAI = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule stats_bam_males:
    input:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        stats = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XY.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"

rule feautreCounts_gene_males:
    input:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        counts = output_directory+"feature_counts_rna/sccREF/{male_sample}_HISAT_geneCounts_XY.txt"
    params:
        featureCounts = featureCounts_path,
        GFF = T2Tv2_annotation_gff_path,
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 2 -t exon -g gene_name -a {params.GFF} -o {output.counts} {input.BAM}"

rule feautreCounts_transcripts_males:
    input:
        BAM = output_directory+"processed_bams/rna/{male_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XY.bam"
    output:
        counts = output_directory+"feature_counts_rna/sccREF/{male_sample}_HISAT_transcriptCounts_XY.txt"
    params:
        featureCounts = featureCounts_path,
        GFF = T2Tv2_annotation_gff_path,
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 2 -t exon -g transcript_id -O -a {params.GFF} -o {output.counts} {input.BAM}"
        
# females
rule samtools_view_females:
    input:
        SAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_XX.sam"
    output:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_XX.bam"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} view -b {input.SAM} > {output.BAM}"

rule bam_sort_females:  
    input:
        IN_BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_XX.bam"
    output:
        sort_BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_XX.bam"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_females:
    input:
        sort_BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_XX.bam"
    output:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam",
        metrics = output_directory+"stats_rna/{female_sample}.XY.picard_mkdup_metrics.txt"
    params:
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} "
        "M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_females:
    input:
        Read_BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_XX.bam"
    output:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam",
    params:
        id = lambda wildcards: config[wildcards.female_sample]["ID"],
        sm = lambda wildcards: config[wildcards.female_sample]["SM"],
        lb = lambda wildcards: config[wildcards.female_sample]["LB"],
        pu = lambda wildcards: config[wildcards.female_sample]["PU"],
        pl = lambda wildcards: config[wildcards.female_sample]["PL"],
        picard = picard_path
    shell:
        "{params.picard} -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} "
        "RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_females:
    input:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        BAI = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam.bai"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} index -in {input.BAM}"

rule stats_bam_females:
    input:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        stats = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_stats_XX.txt"
    params:
        bamtools = bamtools_path
    shell:
        "{params.bamtools} stats -in {input.BAM} > {output.stats}"
 
rule featureCounts_gene_females:
    input:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        counts = output_directory+"feature_counts_rna/sccREF/{female_sample}_HISAT_geneCounts_XX.txt"
    params:
        featureCounts = featureCounts_path,
        GFF = T2Tv2_annotation_gff_path
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 2 -t exon -g gene_name -a {params.GFF} -o {output.counts} {input.BAM}"
        
rule featureCounts_transcripts_females:
    input:
        BAM = output_directory+"processed_bams/rna/{female_sample}_HISAT_pair_trim_sort_mkdup_rdgrp_XX.bam"
    output:
        counts = output_directory+"feature_counts_rna/sccREF/{female_sample}_HISAT_transcriptCounts_XX.txt"
    params:
        featureCounts = featureCounts_path,
        GFF = T2Tv2_annotation_gff_path,
    shell:
        "{params.featureCounts} -T 8 --primary -p -s 2 -t exon -g transcript_id -O -a {params.GFF} -o {output.counts} {input.BAM}"

