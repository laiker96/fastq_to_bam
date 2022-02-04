#! /usr/bin/env python3

def usage():
    '''
    Prints usage message
    '''

    usage_message = '''
--------------------------------------------------------------------------------------------------------|
This script integrates a series of modules for QC, aligning and filtering of fastq.gz files. It expects |
fastq.gz files of the form PREFIX_2.fastq.gz and PREFIX_2.fastq.gz for paired-end experiments and       |
PREFIX.fastq.gz for single-end experiments. PREFIX must NOT include ".".                                |
                                                                                                        |
Directories must be specified without "~" or ".." for referring to home (you may use $HOME) or parent   |
directory. This may cause troubles or end with exceptions or directly not run!).                        |
Use wrapper.py -h/--help for a full description of the program's arguments.                             |
--------------------------------------------------------------------------------------------------------|
'''
    
    print(usage_message)
