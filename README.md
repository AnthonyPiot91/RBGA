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
Some packages are installed in the same environment to simplify the loading of conda environements 
within sbatch scripts.
Others such as chopper, canu and flye need to be installed in separate environments because of 
dependency conflicts. 
Note that creation of the canu environment may take a long time.

```
conda create -n RBGA -c bioconda fastqc trimmomatic a5-miseq nanoplot nanoqc bwa \
polypolish assembly-stats
conda create -n chopper -c bioconda -c conda-forge chopper
conda create -n canu -c bioconda -c conda-forge canu=2.2
conda create -n flye -c bioconda -c conda-forge flye=2.9.2
```

## Install guppy
```

```

## Install snakemake via mamba:
```
mamba install -n HBAW -c bioconda snakemake
```
