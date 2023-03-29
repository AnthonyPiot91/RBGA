# Repetitive Bacterial Genome Assembler (RBGA)

## Description
Pipeline to assemble bacterial genomes with numerous and large repetitive regions using Oxford 
Nanopore long reads and polishing with Illumina Miseq short reads

## Requirements
- The Slurm workload manager
- The conda package manager
- The mamba package manager

## Download this repository
```
git clone https://github.com/AnthonyPiot91/RBGA
```

## Create a conda environment for the workflow and install the required softwares via conda:
```
conda create -n HBAW_1 -c bioconda fastqc trimmomatic a5-miseq nanoplot nanoqc bwa \
polypolish assembly-stats
conda create -n HBAW_2 -c bioconda -c conda-forge chopper
conda create -n canu -c bioconda canu
conda create -n flye -c bioconda flye
 ```

## Install guppy
```

```

## Install snakemake via mamba:
```
mamba install -n HBAW -c bioconda snakemake
```
