[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# Fastq to Bam

This set of scripts conforms a pipeline for processing fastq data from NGS experiments (ATAC-seq, DNase1-seq and ChIP-seq). It performs all processing steps such as adapter removal and mapping to a reference genome. All scripts can be called like independent modules to perform a specific step in the pipeline (e.g., only mapping reads to generate bam file without the need to do adapter-removal previously)


- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [How to run](#how-to-run)

## Overview

DNase-seq counts enriched regions can be detected with the script call_peaks.sh in parallel using the specified number of threads with the MACS2 peak caller. The script autodetects if the BAM files are from a paired-end experiment or a single-end (or mixed parity) experiment and adjust the MACS2 parameters accordingly. Then it estimates the Signal-To-Noise ratio (SNR) of each alignment by calculating the proportion of reads that falls inside peaks (similar to the FrIP score in ATAC-seq experiments). The H3K27ac peaks for each tissue were downloaded from the ENCODE portal as narrowPeak files.

![alt text](pipelineS2.png "pipeline")


## System requirements
### OS requirements
The code has been tested on Ubuntu 20.04.3 LTS

### Dependencies
The following packages are needed to run the code

```
samtools
bowtie2
bwa
bbduk.sh (from the BBMap suite)
Picard
```

### python3 packages dependencies
The following python3 packages are needed to run the code

```
os
```

## Installation
Install from git
```
git clone https://github.com/laiker96/fastq_to_bam
cd fastq_to_bam
```

## How to run
