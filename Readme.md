# General Information

This repository contains code written to analyze brain tumor biopsies.  

# Worflows included

The workflows provided for this omics analysis are created as reproducible Snakemake workflows.  Each directory contains a Snakemake workflow (.snakefile) and a JSON file indicating the names and locations of the data for each sample.  

# Docker for required software

All necessary tools and sex-chromosome complement specific versions of the human reference are packaged into a Docker container, which can be accessed using the following command:
```
docker pull sbplaisier/omics:1.3
```

