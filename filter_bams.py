#! /usr/bin/env python3
import os
import subprocess

# Global variable to specify path to picard MarkDuplicates jar file
picard_exe = ['java', '-jar'
        , '/home/luke/ngsBinaries/picard.jar']

def mark_duplicates(bam_file, threads):
    '''
    Mark duplicates of bam file (optical and PCR are marked) with the 
    corresponding SAM flags bits. It overwrites the original bam file.
    '''

    file_prefix = os.path.splitext(bam_file)[0]

    arguments = [f'I={bam_file}', f'O={bam_file}.flagged'
            , f'M={file_prefix}.metrics'
            , 'ASSUME_SORTED=TRUE']
    subprocess.run(picard_exe + ['MarkDuplicates'] 
            + arguments)
    
    subprocess.run(['mv', f'{bam_file}.flagged', f'{bam_file}'])
    
    arguments_index = ['samtools', 'index'
            , '-@', f'{threads}', f'{bam_file}']
    subprocess.run(arguments_index)


def filter_bam_file(bam_file
        , removedup
        , chromosomes
        , threads
        , MAPQ_threshold
        , paired
        ):
    '''
    Filter bam file with sam flags. If removedupis True, duplicated reads will
    be removed (sam flag bit 1024). Non mapped reads (flag 4) will be removed.
    Non primary alignment are also removed (flag 256)
    '''


    file_prefix = os.path.splitext(bam_file)[0]
    
    if removedup:
        # Bit to remove unmapped, not primary alignment, PCR duplicates
        # and supplementary (chimeric) alignment
        sam_bit_flag_remove = 3332
    else:
        # Bit to remove unmapped, not primary alignment
        # and supplementary (chimeric) alignment
        sam_bit_flag_remove = 2308

    if paired:

        arguments = ['samtools', 'view', '-hb'
                , f'{bam_file}'
                , '--threads', f'{threads}'
                , '-q', f'{MAPQ_threshold}'
                , '-F', f'{sam_bit_flag_remove}'
                , '-f', '3' # Read is paired and mapped in proper paired
                , '-o', f'{file_prefix}_out.bam']
        
        # Fragment size distribution analysis and histogram (.pdf)
        arguments_sizeMetrics = ['-I', f'{bam_file}'
                , '-O', f'{file_prefix}_sizes.txt'
                , '-H', f'{file_prefix}_sizes.pdf']
        subprocess.run(picard_exe + ['CollectInsertSizeMetrics'] 
            + arguments_sizeMetrics)

    else:

        arguments = ['samtools', 'view', '-hb'
                , f'{bam_file}'
                , '--threads', f'{threads}'
                , '-q', f'{MAPQ_threshold}'
                , '-F', f'{sam_bit_flag_remove}'
                , '-o', f'{file_prefix}_out.bam']

    #chromosomes = chromosomes.split(",")
    subprocess.run(arguments + chromosomes)
    
    arguments_index = ['samtools', 'index'
            , '-@', f'{threads}', f'{file_prefix}_out.bam']
    subprocess.run(arguments_index)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    #Add arguments to parser
    parser.add_argument("-b", "--bam", 
                    help = '''bam file to process''',
                    type = str)

    parser.add_argument("-r", "--remove-dups",
                    help = '''bool to remove duplicate reads''',
                    type = str,
                    default = False)

    parser.add_argument("-c", "--chromosomes",
                    help = '''chomosomes to keep''',
                    type = str,
                    default = None)

    parser.add_argument("-t", "--threads",
                    help = '''Number of threads to use''',
                    type = str,
                    default = 1)
                    
    parser.add_argument("-m", "--mapq",
                    help = '''MAPQ threshold to use''',
                    type = int,
                    default = 20)
                    
    parser.add_argument("-p", "--paired",
                    help = '''Files are paired or unpaired''',
                    type = bool,
                    default = False)



    args = parser.parse_args()
    args = (vars(args))
    args = list(args.values())
    
    args[2] = args[2].split(",")
    mark_duplicates(args[0], args[3])

    filter_bam_file(*args)
