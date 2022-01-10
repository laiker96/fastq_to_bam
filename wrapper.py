#! /usr/bin/env python3

import os
import sys
import json
import argparse
import subprocess
import usage
import fastq_2_bam
import create_jsons


#Create argument parser for main pipeline
parser = argparse.ArgumentParser(usage = usage.usage())
#Add arguments to parser
parser.add_argument("-w", "--working_directory", 
        help = '''Directory with fastq.gz files to process''',
        type = str)


parser.add_argument("-a", "--adapters",
        help = '''adapter fasta file with sequences to scan reads 
        for adapter trimming''',
        type = str)


parser.add_argument("-t", "--threads",
        help = "number of threads to use for ALL analyses",
        type = int,
        default = 1)


parser.add_argument("-m", "--minlength",
        help = '''minimum fraction of original read length for reads to 
        hold after adapter trimming. Reads shorter than this fraction value 
        WILL be discarded''',
        type = float,
        default = 0.5)


parser.add_argument("-i", "--index",
        help = "index directory (PE) or prefix (SE)",
        type = str)


parser.add_argument("-c", "--chromosomes",
        help = '''chromosomes to keep in bam file (comma separated). 
        If not specified, all chromosomes will be kept''',
        type = str)


parser.add_argument("-Q", "--MAPQ_threshold",
        help = "MAPQ threshold value for bam file filtering",
        type = int)


parser.add_argument("-r", "--removedup",
        help = "flag for removing duplicates from bam files",
        action = "store_true")


parser.add_argument("-q", "--quality_dir",
        help = "name of quality reports directory",
        type = str,
        default = 'quality')


parser.add_argument("-A", "--aligner",
        help = "aligner program (bwa mem or bowtie2)",
        type = str,
        default = 'bwa')


parser.add_argument("-o", "--output_dir",
        help = '''name of output directory (trimmed fastq.gz and bam 
        files will be written here)''',
        type = str,
        default = 'output_files')


parser.add_argument("-v", "--verbose",
        help = "Print verbose output",
        action = "store_true")

args = parser.parse_args()

directory = args.working_directory
threads = args.threads
adapters = args.adapters
minlength = args.minlength
index = args.index
chromosomes = (args.chromosomes.split(",") 
        if args.chromosomes is not None else [])
mapq = args.MAPQ_threshold
removedup = args.removedup
qualdir = args.quality_dir
outdir = args.output_dir
aligner = args.aligner
verbose = args.verbose

def summarize_directory(
        directory
        , threads
        , adapters
        , minlength
        , index
        , chromosomes
        , mapq
        , removedup
        , qualdir
        , outdir
        , verbose
        ):
    '''
    Create .json argument files for all fastq files in the given directory    
    '''

    prefixes = create_jsons.scan_directory(directory)
    SE = prefixes[0]
    PE = prefixes[1]
    json_list = []
    args = [threads
            , adapters
            , minlength
            , index
            , chromosomes
            , mapq
            , removedup
            , qualdir
            , outdir
            , verbose]

    for prefix in SE:
        
        json_list.append(create_jsons.create_jsons(prefix, False, *args))

    for prefix in PE:

        json_list.append(create_jsons.create_jsons(prefix, True, *args))
    
    return json_list

def main(directory
            , threads
            , adapters
            , minlength
            , index
            , chromosomes
            , mapq
            , removedup
            , qualdir
            , outdir
            , aligner
            , verbose):
    '''
    Main function
    '''

    json_list = summarize_directory(directory
            , threads
            , adapters
            , minlength
            , index
            , chromosomes
            , mapq
            , removedup
            , qualdir
            , outdir
            , verbose)

    for json_file in json_list:
        argument_dict = fastq_2_bam.read_json(json_file)
        fastq_2_bam.main(argument_dict['qualdir']
            , argument_dict['outdir']
            , argument_dict['files']
            , argument_dict['threads']
            , argument_dict['paired']
            , argument_dict['adapters']
            , argument_dict['minlength']
            , argument_dict['index']
            , argument_dict['removedup']
            , argument_dict['mapq']
            , argument_dict['chromosomes']
            , aligner
            , argument_dict['verbose']
            )

main(directory
            , threads
            , adapters
            , minlength
            , index
            , chromosomes
            , mapq
            , removedup
            , qualdir
            , outdir
            , aligner
            , verbose)


