# Repository Contents

This repository contains code written to analyze brain tumor biopsies that were processed using RNA sequencing for gene expression measurements and whole exome sequencing (WES) for DNA variant calling. Each directory within this repository includes a readme that details how the workflows should be run. 

# Worflows included

The workflows provided for this omics analysis are created as reproducible Snakemake workflows.  Each directory contains a Snakemake workflow (.snakefile) and a JSON file indicating the names and locations of the data for each sample.  

# Docker for required software

All necessary tools and sex-chromosome complement specific versions of the human reference are packaged into a Docker container, which can be accessed using the following command:
```
docker pull sbplaisier/omics:1.3
```

