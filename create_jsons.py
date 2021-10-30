import json
import os
import sys

def scan_directory(directory):
    '''
    Scan directory and returns tuple of lists with
    single-end and paired-end prefixes for creating
    .json files
    '''
    
    # Create list of files paths
    if directory == '.':
        files = os.listdir()
        files_path = files
    else:
        files = os.listdir(directory)
        files_path = [os.path.join(directory, file) for file in files]
    # Filter for .fastq.gz files
    fastq_files = [file for file in files_path 
            if file.endswith(".fastq.gz")]
    
    # Split in PE and SE files

    paired_end_files = [fastq for fastq in fastq_files 
            if ('_1.fastq.gz' in os.path.basename(fastq) 
                or '_2.fastq.gz' in os.path.basename(fastq))]

    single_end_files = [fastq for fastq in fastq_files 
            if fastq not in paired_end_files] 
    
    # Create prefix files lists
    single_end_prefix = ['.'.join(fastq.split('.')[:-2])
            for fastq in single_end_files]

    paired_end_prefix = ['_'.join(fastq.split('_')[:-1])
            for fastq in paired_end_files]
    paired_end_prefix = list(set(paired_end_prefix))
    
    # Check that all PE files have their respective pair file
    checker_PE = [((prefix + '_1.fastq.gz') in paired_end_files 
        and (prefix + '_2.fastq.gz') in paired_end_files)
            for prefix in paired_end_prefix]
    
    # Exit program if there are inconsistencies between file names of
    # PE files
    if not all(checker_PE):
        print('Inconsistent names of paired-end files!')
        sys.exit(1)
    

    return single_end_prefix, paired_end_prefix

def create_jsons(
        file_prefix
        , paired
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
        , json_output_name = None
        ):
    '''
    Create .json files for processing data. The files are writen
    to disk
    '''
    
    if json_output_name is None:
        json_output_name = f'{file_prefix}_values.json'
    
    # Initialize output file
    dict_values = dict()
    
    if paired:
        files = [file_prefix + '_1.fastq.gz', file_prefix + '_2.fastq.gz']
    else:
        files = [file_prefix + '.fastq.gz']
    
    # Keys and values definition
    dict_values['files'] = files
    dict_values['paired'] = paired
    dict_values['threads'] = threads
    dict_values['adapters'] = adapters
    dict_values['minlength'] = minlength
    dict_values['index'] = index
    dict_values['chromosomes'] = chromosomes
    dict_values['mapq'] = mapq
    dict_values['removedup'] = removedup
    dict_values['qualdir'] = qualdir
    dict_values['outdir'] = outdir
    dict_values['verbose'] = verbose
    
    # Write dict object to disk (serialize to .json format)
    with open(json_output_name, 'wt') as output:

        json_file = json.dumps(dict_values, indent = 4)
        output.write(json_file)

    # Return name of file for subsequent analyses
    return json_output_name


