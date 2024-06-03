import os

configfile: "Trimmed_FastqFilePaths_round2.rdgroups.withnormals.json"

strelka_path = "/home/splaisie/software/Strelka/strelka-2.9.10.centos6_x86_64/bin"
gatk_path = "gatk"

refpath_XX ="T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
refpath_XY = "T2Tv2_sex_complement_reference/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",

rule all:
   input: # variant calling and filtering with strelka
        expand("variant_calling_strelka/{sample}/runWorkflow.py", sample=config["female_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XX.snvs.vcf.gz", sample=config["female_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XX.snvs.pass.vcf.gz", sample=config["female_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XX.indels.vcf.gz", sample=config["female_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XX.indels.pass.vcf.gz", sample=config["female_tumors"]),
        expand("variant_calling_strelka/{sample}/runWorkflow.py", sample=config["male_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XY.snvs.vcf.gz", sample=config["male_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XY.snvs.pass.vcf.gz", sample=config["male_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XY.indels.vcf.gz", sample=config["male_tumors"]),
        expand("variant_calling_strelka/{sample}/results/variants/somatic.XY.indels.pass.vcf.gz", sample=config["male_tumors"])

# call females
rule configure_tumor_with_matched_normal_XX:
    input:
        ref = refpath_XX,
        normal_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XX.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XX.dedup.bam")
    output:
        "variant_calling_strelka/{sample}/runWorkflow.py"
    params:
        strelka = strelka_path,
        run_dir = "variant_calling_strelka/{sample}",
    shell:
        """
        {params.strelka}/configureStrelkaSomaticWorkflow.py --referenceFasta {input.ref} --tumorBam {input.tumor_bam} --normalBam {input.normal_bam} --runDir {params.run_dir} --exome
        """

rule run_tumor_with_matched_normal_XX:
    input:
        ref = refpath_XX,
        normal_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XX.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XX.dedup.bam"),
        workflow = "variant_calling_strelka/{sample}/runWorkflow.py"
    output:
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.XX.snvs.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.XX.indels.vcf.gz"
    params:
        run = "variant_calling_strelka/{sample}/runWorkflow.py",
    shell:
        """
        {params.run} -m local -j 20
        """

rule select_pass_variants_XX:
    input:
        ref = refpath_XX,
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.XX.snvs.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.XX.indels.vcf.gz"
    output:
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.XX.snvs.pass.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.XX.indels.pass.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.snvs} --exclude-filtered -O {output.snvs}
        {params.gatk} SelectVariants -R {input.ref} -V {input.indels} --exclude-filtered -O {output.indels}
        """

# call males
rule configure_tumor_with_matched_normal_XY:
    input:
        ref = refpath_XY,
        normal_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.dedup.bam")
    output:
        "variant_calling_strelka/{sample}/runWorkflow.py"
    params:
        strelka = strelka_path,
        run_dir = "variant_calling_strelka/{sample}",
    shell:
        """
        {params.strelka}/configureStrelkaSomaticWorkflow.py --referenceFasta {input.ref} --tumorBam {input.tumor_bam} --normalBam {input.normal_bam} --runDir {params.run_dir} --exome
        """

rule run_tumor_with_matched_normal_XY:
    input:
        ref = refpath_XY,
        normal_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["normal"] + ".T2Tv2.XY.dedup.bam"),
        tumor_bam = lambda wildcards: os.path.join("alignments/", config[wildcards.sample]["ID"] + ".T2Tv2.XY.dedup.bam"),
        workflow = "variant_calling_strelka/{sample}/runWorkflow.py"
    output:
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.XY.snvs.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.XY.indels.vcf.gz"
    params:
        run = "variant_calling_strelka/{sample}/runWorkflow.py",
    shell:
        """
        {params.run} -m local -j 20
        """

rule select_pass_variants_XY:
    input:
        ref = refpath_XY,
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.XY.snvs.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.XY.indels.vcf.gz"
    output:
        snvs = "variant_calling_strelka/{sample}/results/variants/somatic.XY.snvs.pass.vcf.gz",
        indels = "variant_calling_strelka/{sample}/results/variants/somatic.XY.indels.pass.vcf.gz"
    params:
        gatk = gatk_path
    shell:
        """
        {params.gatk} SelectVariants -R {input.ref} -V {input.snvs} --exclude-filtered -O {output.snvs}
        {params.gatk} SelectVariants -R {input.ref} -V {input.indels} --exclude-filtered -O {output.indels}
        """
