# fastq_to_bam

This set of scripts conforms a pipeline for processing fastq data from NGS experiments (ATAC-seq, DNase1-seq and ChIP-seq). It performs all processing steps such as adapter removal and mapping to a reference genome. All scripts can be called like independent modules to perform a specific step in the pipeline (e.g., only mapping reads to generate bam file without the need to do adapter-removal previously)


