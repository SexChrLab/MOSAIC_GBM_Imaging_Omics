import os

configfile: "Trimmed_FastqFilePaths_round2.rdgroups.withnormals.json"

gatk_path = "gatk"

refpath_XX = "/references/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
refpath_XY = "/references/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",

working_dir = "/data/Test_variant_callers/"
align_dir = "/data/Test_variant_callers/alignments/"

rule all:
   input: # variant calling and filtering with MuTect2
        expand(working_dir+"variant_calling_mutect2/{sample}.somatic.vcf.gz", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_mutect2/{sample}.somatic.filtered.vcf.gz", sample=config["female_tumors"]),
        expand(working_dir+"variant_calling_mutect2/{sample}.somatic.filtered.pass.vcf.gz", sample=config["female_tumors"]),
		expand(working_dir+"variant_calling_mutect2/{sample}.somatic.vcf.gz", sample=config["male_tumors"]),
        expand(working_dir+"variant_calling_mutect2/{sample}.somatic.filtered.vcf.gz", sample=config["male_tumors"]),
        expand(working_dir+"variant_calling_mutect2/{sample}.somatic.filtered.pass.vcf.gz", sample=config["male_tumors"])

rule tumor_with_matched_normal_XY:
    input:
        ref = refpath_XY,
        normal_bam = lambda wildcards: os.path.join(align_dir, config[wildcards.sample]["normal"] + ".T2Tv2.XY.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join(align_dir, config[wildcards.sample]["ID"] + ".T2Tv2.XY.dedup.bam")
    output:
        os.path.join(working_dir+"variant_calling_mutect2/", "{sample}.somatic.vcf.gz")
    params:
        gatk = gatk_path,
        sm = lambda wildcards: config[wildcards.sample]["normal"]
    shell:
        """
        {params.gatk} Mutect2 -R {input.ref} -I {input.tumor_bam} -I {input.normal_bam} -normal {params.sm} -O {output}
        """

rule filter_XY:
    input:
        ref = refpath_XY,
        unfiltered = os.path.join(working_dir+"variant_calling_mutect2/", "{sample}.somatic.vcf.gz")
    output:
        filtered = os.path.join(working_dir+"variant_calling_mutect2/", "{sample}.somatic.filtered.vcf.gz")
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} FilterMutectCalls -R {input.ref} -V {input.unfiltered} -O {output.filtered}
        """

rule select_pass_variants_XY:
    input:
        ref = refpath_XY,
        vcf = os.path.join(working_dir+"variant_calling_mutect2/", "{sample}.somatic.filtered.vcf.gz")
    output:
        os.path.join(working_dir+"variant_calling_mutect2/", "{sample}.somatic.filtered.pass.vcf.gz")
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.vcf} --exclude-filtered -O {output}
        """
