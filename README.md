# Repetitive Bacterial Genome Assembler (RBGA)

## Description
Pipeline to assemble bacterial genomes with numerous and large repetitive regions using Oxford 
Nanopore long reads and polishing with Illumina Miseq short reads

## Requirements
- The Slurm workload manager
- The conda package manager

## Download this repository
```
git clone https://github.com/AnthonyPiot91/RBGA
```

## Install the required softwares via conda:
Some packages are installed in the same environment to simplify the loading of conda environements 
within sbatch scripts.

Others such as chopper and flye need to be installed in separate environments because of 
dependency conflicts. 

```
conda create -n RBGA \
	     -c bioconda \
	     -c conda-forge \
	     fastqc \
	     trimmomatic \
	     nanoplot \
	     nanoqc \
	     bwa \
	     polypolish \
	     assembly-stats \
	     seqtk
conda create -n chopper \
	     -c bioconda \
	     -c conda-forge \
	     chopper
conda create -n flye \
	     -c bioconda \
	     -c conda-forge \
	     flye=2.9.2
```
