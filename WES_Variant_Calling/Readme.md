# Variant Calling Workflow

This set of workflows are given to perform somatic variant calling using three algorithms.  Each algorithm strikes a different balance between sensitivity and specificity so we intend to take the intersection between multiple algorithms to find identify variants.  

# Files included

Trimmed_FastqFilePaths_round2.rdgroups.withnormals.json: most recent config json containing sample names of all samples to run after trimming

variant_calling_MuTect2.snakefile: workflow for variant calling using the MuTect2 package

variant_calling_Varscan.snakefile: workflow for variant calling using the Varscan package

variant_calling_strelka.snakefile: workflow for variant calling using the Strelka package

# Running workflow

For dry run: 
```
snakemake -np -s variant_calling_[MuTect2|Varscan|strelka].snakefile —-latency-wait=60
```

For full execution: 
```
snakemake -s variant_calling_[MuTect2|Varscan|strelka].snakefile —-latency-wait=60
```

You can include the -c4 Snakemake flag for parallelization of jobs if your system supports this.  
