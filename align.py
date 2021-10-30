#! /usr/bin/env python3
import subprocess
import os
import shutil


def align_reads(files, outdir, paired, index, threads, aligner):
    '''
    Align reads to reference genome using index files prefix (SE) or 
    directory(PE) and a list of fastq.gz files names. 
    '''
    print(files)
    print(paired)
    if paired:
        # Save file prefix to variable (files must be .fastq.gz)
        file_prefix = '.'.join(files[0].split('.')[:-2])
        # Remove the underscore from the prefix
        file_prefix = file_prefix[:-2] if file_prefix[-2] == "_" else file_prefix
        # New prefix name (for sam and bam files)
        new_name = os.path.join(outdir, os.path.basename(file_prefix))
        
        # Arguments definition for bwa aligner
        if aligner == "bwa":
            arguments = ['bwa', 'mem', '-t'
                    , f'{threads}', f'{index}'
                    , '-o', f'{file_prefix}.sam'
                    , f'{files[0]}', f'{files[1]}']
        # Call aligner
        else:
            arguments = ['bowtie2', '--very-sensitive', '--threads'
                    , f'{threads}', '-x', f'{index}'
                    , '-S', f'{file_prefix}.sam'
                    , '-1', f'{files[0]}', '-2',f'{files[1]}']
            
        subprocess.run(arguments)

    else:
        # Save file prefix to variable (files must be .fastq.gz)
        file_prefix = '.'.join(files[0].split('.')[:-2])
        # New prefix name (for sam and bam files)
        new_name = os.path.join(outdir, os.path.basename(file_prefix))
        
        # Arguments definition for bwa aligner
        if aligner == "bwa":
            arguments = ['bwa', 'mem', '-t'
                    , f'{threads}', f'{index}'
                    , '-o', f'{file_prefix}.sam'
                    , f'{files[0]}']
        else:
            arguments = ['bowtie2', '--very-sensitive', '--threads'
                    , f'{threads}', '-x', f'{index}'
                    , '-S', f'{file_prefix}.sam'
                    , '-U', f'{files[0]}']
            
        # Call aligner
        subprocess.run(arguments)
        
    # Arguments definition for samtools (convert sam to bam)
    arguments_samtools = ['samtools', 'view', '-hb'
            , '--threads', f'{threads}'
            , '-o', f'{new_name}.bam'
            , f'{file_prefix}.sam']
    # Call samtools
    subprocess.run(arguments_samtools)
    
    # Remove sam file (keep only bam)
    os.remove(f'{file_prefix}.sam')
    
    # Return bam file name for subsequent analyses
    return f'{new_name}.bam'

def sort_and_index(bam_file, threads):
    '''
    Sort and index bam file with samtools
    '''

    # Arguments definition for samtools sort
    arguments = ['samtools', 'sort', '--threads'
            , f'{threads}', f'{bam_file}'
            , '-o', f'{bam_file}.sorted']
    # Call samtools sort
    subprocess.run(arguments)
    
    # Overwrite unsorted bam with sorted bam file
    subprocess.run(['mv', f'{bam_file}.sorted', f'{bam_file}'])
    
    # Arguments definition for samtools index
    arguments_index = ['samtools', 'index'
            , '-@', f'{threads}', f'{bam_file}']
    # Call samtools index
    subprocess.run(arguments_index)
