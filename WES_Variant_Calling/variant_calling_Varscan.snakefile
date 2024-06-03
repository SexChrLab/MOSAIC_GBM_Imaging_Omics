import os

configfile: "Trimmed_FastqFilePaths_round2.rdgroups.withnormals.json"

varscan_path = "/home/splaisie/software/Varscan/VarScan.v2.3.9.jar"
gatk_path = "gatk"

working_dir = "/home/ext_hawkins_daarud_andrea_mayo_e/"
align_dir = "/home/ext_hawkins_daarud_andrea_mayo_e/FuseMount/hawkinsdaarud/WES/alignments/"

refpath_XX = "/references/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
refpath_XY = "/references/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",

rule all:
   input: # variant calling and filtering with varscan
        expand(working_dir+"variant_calling_varscan/pileups/{sample}.T2Tv2.XX.pileup", sample=config["female_samples"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.XX.snp", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.snp.Somatic.XX.hc", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.snp.Somatic.XX.hc.filter", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.snp.Somatic.XX.hc.filter.vcf", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.XX.indel", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_varscan/pileups/{sample}.T2Tv2.XY.pileup", sample=config["male_samples"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.XY.snp", sample=config["male_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.snp.Somatic.XY.hc", sample=config["male_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.snp.Somatic.XY.hc.filter", sample=config["male_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.snp.Somatic.XY.hc.filter.vcf", sample=config["male_tumors"]),
        expand(working_dir+"variant_calling_varscan/{sample}.varscan.XY.indel", sample=config["male_tumors"])

rule bam_pileup_XX: #for both normal and tumor
    input:
        ref = refpath_XX,
        bam = os.path.join(alignment_directory,"alignments/", "{sample}.T2Tv2.XX.dedup.bam")
    output:
        temp("variant_calling_varscan/pileups/{sample}.T2Tv2.XX.pileup")
    threads: 4
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} > {output}
        """

rule run_varscan_XX:
    input:
        ref = refpath_XX,
        normal_pileup = lambda wildcards: os.path.join("variant_calling_varscan/pileups/", config[wildcards.sample]["normal"] + ".T2Tv2.XX.pileup"),
        tumor_pileup = lambda wildcards: os.path.join("variant_calling_varscan/pileups/", config[wildcards.sample]["ID"] + ".T2Tv2.XX.pileup"),
    output:
        snp = "variant_calling_varscan/{sample}.varscan.snp",
        indel = "variant_calling_varscan/{sample}.varscan.indel"
    params:
        varscan = varscan_path,
        basename = "variant_calling_varscan/{sample}.varscan"
    shell:
        """
        java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05
        """

rule bam_pileup_XY: #for both normal and tumor
    input:
        ref = refpath_XY,
        bam = os.path.join(alignment_directory,"alignments/", "{sample}.T2Tv2.XY.dedup.bam")
    output:
        temp("variant_calling_varscan/pileups/{sample}.T2Tv2.XY.pileup")
    threads: 4
    shell:
        """
        samtools mpileup -f {input.ref} {input.bam} > {output}
        """

rule run_varscan_XY:
    input:
        ref = refpath_XY,
        normal_pileup = lambda wildcards: os.path.join("variant_calling_varscan/pileups/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.pileup"),
        tumor_pileup = lambda wildcards: os.path.join("variant_calling_varscan/pileups/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.pileup"),
    output:
        snp = "variant_calling_varscan/{sample}.varscan.snp",
        indel = "variant_calling_varscan/{sample}.varscan.indel"
    params:
        varscan = varscan_path,
        basename = "variant_calling_varscan/{sample}.varscan"
    shell:
        """
        java -jar {params.varscan} somatic {input.normal_pileup} {input.tumor_pileup} {params.basename} –min-coverage 10 –min-var-freq 0.08 –somatic-p-value 0.05
        """

rule isolate_calls_by_type_and_confidence:
    input: 
        snp = "variant_calling_varscan/{sample}.varscan.snp",
    output: 
        snp_somatic_hc = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc",
    params: 
        varscan = varscan_path
    shell:
        """
        java -jar {params.varscan} processSomatic {input.snp}
        """

rule somatic_filter:
    input: 
        snp_somatic_hc = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc",
        indel = "variant_calling_varscan/{sample}.varscan.indel",
    output: 
        snp_somatic_hc_filter = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter",
    params: 
        varscan = varscan_path
    shell:
        """
        java -jar {params.varscan} somaticFilter {input.snp_somatic_hc} -indel-file {input.indel} -output-file {output.snp_somatic_hc_filter}
        """

rule convert_to_vcf:
    input:
        snp_somatic_hc_filter = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter",
    output:
        snp_somatic_hc_filter_vcf = "variant_calling_varscan/{sample}.varscan.snp.Somatic.hc.filter.vcf",
    shell:
        """
        cat {input.snp_somatic_hc_filter} | awk '{{print $1"\t" $2"\t" "." "\t" $3 "\t" $4}}' > {output.snp_somatic_hc_filter_vcf}
        """

