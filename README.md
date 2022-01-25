[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# Fastq to Bam

This set of scripts conforms a pipeline for processing fastq data from NGS experiments (ATAC-seq, DNaseI-seq and ChIP-seq). It performs all processing steps such as adapter removal and mapping to a reference genome. All scripts can be called like independent modules to perform a specific step in the pipeline (e.g., only mapping reads to generate bam file without the need to do adapter-removal previously)


- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [How to run](#how-to-run)

## Overview


```preprocess_fastq.py``` runs Fastqc on the fastq.gz files and bbduk.sh with a fasta of adapters for adapter removal and trimming. It then runs Fastqc on the trimmed fastq.gz files an save the reports as .html files. You should check this files.

```align.py``` runs bwa or bowtie2 aligner with the given reference genome index and the trimmed fastq.gz files. It also compress the sam output file to a bam file and indexes it.

```filter_bam.py``` filters the bam files with MAPQ value and sam flag bits.

The ```wrapper.py``` script automatically processes all files in the given working directory, running all steps of the pipeline. All fastq files must be named as ```<fileID>.fastq.gz``` for SE files and ```<fileID>_1.fastq.gz <fileID>_1.fastq.gz``` for PE files.


## System requirements
### OS requirements
The code has been tested on Ubuntu 20.04.3 LTS with python version 3.8.8

### Dependencies
The following packages are needed to run the code

```
samtools
bowtie2
bwa
bbduk.sh (from the BBMap suite)
Picard
```


## Installation
Install from git
```
git clone https://github.com/laiker96/fastq_to_bam
cd fastq_to_bam
```

## How to run
To run each module independently:

Preprocess fastq files for alignment
```preprocess_fastq.py -f <FASTQ.GZ_FILES> -t <THREADS> -o <OUTDIR_NAME> -p <PAIRED_BOOL> -a <ADAPTER_FASTA_FILE> -m <MINIMUM_FRACTION_LENGTH>```
