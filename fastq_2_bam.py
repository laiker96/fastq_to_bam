#! /usr/bin/env python3

import os
import sys
import json
import subprocess
import usage
import preprocess_fastq
import align
import filter_bams



def read_json(inputJSON):
    '''
    Read json files from disk and save them to dict objects
    '''

    with open(inputJSON, 'rt') as json_file:

        argument_dict = json.load(json_file)
    
    return argument_dict


# Create output directory if not in path
def create_dirs(quality_dir, output_dir, verbose):
    '''
    Create directories for quality reports and output files if these are not 
    in path. Returns a tuple with the names of these directories
    '''

    if not os.path.isdir(quality_dir):
        os.mkdir(quality_dir)
        if verbose:
            print('Quality data directory created')
    else:
        if verbose:
            print('Quality directory already exists!')
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        if verbose:
            print('Output data directory created')
    else:
        if verbose:
            print('Output directory already exists!')

    return quality_dir, output_dir

def main(qualdir
        , outdir
        , files
        , threads
        , paired
        , adapters
        , minlength
        , index
        , removedup
        , mapq
        , chromosomes
        , aligner
        , verbose
        ):
    '''
    Main pipeline function
    '''
    # Create directories if not in path
    create_dirs(qualdir, outdir, verbose)
    # QC analysis of fastq files
    preprocess_fastq.fastq_qc(files, threads, qualdir)
    # Adapter trimming
    archivos = preprocess_fastq.trim_adapters(files, threads, outdir
            , paired, adapters, minlength)
    # QC analysis of processed fastq files
    preprocess_fastq.fastq_qc(archivos, threads, qualdir)
    # Align reads to reference with bowtie2 (SE) or STAR (PE)
    print(archivos)
    out = align.align_reads(archivos, outdir, paired, index, threads, aligner)
    # Sort and index bam files
    align.sort_and_index(out, threads)
    # Mark duplicated reads with picard MarkDuplicates
    filter_bams.mark_duplicates(out, threads)
    # Filter bam files
    filter_bams.filter_bam_file(out, removedup, chromosomes, threads, mapq, paired)
