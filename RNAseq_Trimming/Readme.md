# Trimming Workflow

To remove low quality RNA sequencing reads, we used trimming software bbduk from the bbmap software package to trim RNA sequencing data files.

# Files included

FastqFilePaths.json: most recent config json containing sample names of all samples to run (pulling out the sample where FastQC indicated that the file was truncated)

trim.snakemake: snakemake file that reads the fastq file path, forward reads file, and reverse reads file of samples listed in config json and runs the bbduk command with trimming parameters

# Running workflow

For dry run: 
```
snakemake -np -s trim.snakefile —-latency-wait=60
```

For full execution: 
```
snakemake -s trim.snakefile —-latency-wait=60
```

You can include the -c4 Snakemake flag for parallelization of jobs if your system supports this.  Upon completion, trimmed fastq files will be in trimmed_fastqs/ directory.
