# RNA Sequencing Workflow

This workflow includes sex-informed alignment using the HISAT2 algorithm to the human telomere-to-telomere reference genome (CHM13v2) and quantification of gene expression using featureCounts.

# Files included

FastqFilePaths.json: most recent config json containing sample names of all samples to run

rnaseq_data_processing.snakemake: snakemake workflow for alignment and gene expression quantification

# Running workflow

For dry run: 
```
snakemake -np -s rnaseq_data_processing.snakefile —-latency-wait=60
```

For full execution: 
```
snakemake -s rnaseq_data_processing.snakefile —-latency-wait=60
```

You can include the -c4 Snakemake flag for parallelization of jobs if your system supports this.  Upon completion, count files will be placed in feature_counts_rna/ directory.
