configfile: "FastqFilePaths_local.json"

# Illumina TruSeq adapters
adapter_path = "adapter_sequence.fa"

# Trim with BBDuk trimming software from the Bbmap package from the Joint Genome Institute (https://sourceforge.net/projects/bbmap/)
#   install bbmap: conda install -c bioconda bbmap
#   test install: enter 'bbduk.sh' on command line in conda environment you installed from, should display all the options for bbduk

rule all:
    input:
        expand("trimmed_fastqs/{sample}_trimmed_{run}.fq.gz", sample=config["sample_names"],run=["1","2"])

rule trim_adapters_paired_bbduk_dna:
    input:
        fq1 = lambda wildcards: os.path.join(config[wildcards.sample]["fq_path"], config[wildcards.sample]["fq1"]),
        fq2 = lambda wildcards: os.path.join(config[wildcards.sample]["fq_path"], config[wildcards.sample]["fq2"])
    output:
        out_fq1 = "trimmed_fastqs/{sample}_trimmed_1.fq.gz",
        out_fq2 = "trimmed_fastqs/{sample}_trimmed_2.fq.gz"
    params:
        adapter = adapter_path
    threads:
        2
    shell:
        "bbduk.sh -Xmx3g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref={params.adapter} qtrim=rl trimq=10 minlen=75 maq=20"

