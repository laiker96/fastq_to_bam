#! /usr/bin/env python3
import os
import subprocess


def get_picard_config_path():
    '''
    Build picard config file
    '''

    script_path = os.path.realpath(__file__)
    script_dir = os.path.dirname(script_path)
    picard_config_path = os.path.join(script_dir, '.picard_path.conf')
    
    if not os.path.exists(picard_config_path):

        picard_path = input('Please enter full path to picard .jar file: ')
        
        with open(picard_config_path, 'wt') as conf_file:
            conf_file.write(picard_path)
    
    with open(picard_config_path, 'rt') as conf_file:
        picard_path = conf_file.readline()
    
        
    return ['java', '-jar', picard_path]

#Global variable to access the picard jar file
picard_exe = get_picard_config_path()


def mark_duplicates(bam_file, threads):
    '''
    Mark duplicates of bam file (optical and PCR are marked) with the 
    corresponding SAM flags bits. It overwrites the original bam file.
    '''

    file_prefix = os.path.splitext(bam_file)[0]

    arguments = [f'I={bam_file}', f'O={bam_file}.flagged'
            , f'M={file_prefix}.metrics'
            , 'ASSUME_SORTED=TRUE', 'MAX_RECORDS_IN_RAM=10000000']
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
    Filter bam file with sam flags. If removedup is True, duplicated reads will
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
                , '-f', '3' # Read is paired and mapped in proper pair
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
                    help = '''Remove duplicate reads from alignment file''',
                    action = "store_true")

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
                    help = '''Alignment files are paired''',
                    action = 'store_true')



    args = parser.parse_args()
    args = (vars(args))
    args = list(args.values())
    
    args[2] = args[2].split(",")
    print(args)
    mark_duplicates(args[0], args[3])

    filter_bam_file(*args)
