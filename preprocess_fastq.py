#! /usr/bin/env python3
import os
import subprocess
import shutil

def fastq_qc(files
        , threads
        , output_dir
        ):
    '''
    Quality check analysis of fastq files with FastQC. Reports (zip and htmls)
    files are saved in output_dir
    '''
    
    # Call FastQC
    subprocess.run(['fastqc', '--threads', f'{threads}'] + files)
    
    # Move QC reports to quality directory
    for file in files:
        
        # Create file prefix (file name without .fastq.gz)
        file_prefix = ''.join(file.split('.')[:-2])
        
        print(file_prefix)
        # Current QC zip and htmls file names
        html_name, zip_name = (file_prefix + '_fastqc.html'
                , file_prefix + '_fastqc.zip')

        # New names based on quality directory
        html_new, zip_new = (os.path.join(output_dir, os.path.basename(html_name))
                , os.path.join(output_dir, os.path.basename(zip_name)))
        
        # Rename files to desired location
        shutil.move(html_name, html_new)
        shutil.move(zip_name, zip_new)


def trim_adapters(files
        , threads
        , outdir
        , paired
        , adapters
        , minlength
        ):
    '''
    Trim adapters from fastq files using bbduk. Output trimmed fastq.gz files
    are saved in outdir directory. It returns a list of the names of the 
    trimmed fastq.gz files for subsequent analysis
    '''
    
    if paired:
        
        # Define basenames for each pair of fastq.gz files
        file_basenames = os.path.basename(files[0]), os.path.basename(files[1])
        # Define global name for experiment
        file_basename = os.path.basename(files[0]).split('_')[0]
        
        # Arguments definition for bbduk.sh
        arguments = ['Xmx12g', f'in1={files[0]}'
                , f'in2={files[1]}'
                , f'out1={outdir}/trim_{file_basenames[0]}'
                , f'out2={outdir}/trim_{file_basenames[1]}'
                , f'ref={adapters}'
                , 'prealloc=t', 'ktrim=r', 'k=23', 'hdist=1'
                , 'tbo=t', 'tpe=t', 'mink=11', f'mlf={minlength}'
                , 'rcomp=t', f'stats={outdir}/{file_basename}.trimlog']
        
        # The funcion returns the names of the trimmed fastq.gz files
        return_value = [f'{outdir}/trim_{file_basenames[0]}', f'{outdir}/trim_{file_basenames[1]}']
        # return_value = [value[2:] if value[:2] == './' else value for value in return_value]
    
    else:
        
        # Define global name for experiment
        file_basename = os.path.basename(files[0])
        
        # Arguments definition for bbduk.sh
        arguments = ['Xmx12g', f'in={files[0]}'
                , f'out={outdir}/trim_{file_basename}'
                , f'ref={adapters}'
                , 'prealloc=t', 'ktrim=r', 'k=23', 'hdist=1'
                , 'mink=11', f'mlf={minlength}', 'rcomp=t'
                , f'stats={outdir}/{os.path.splitext(file_basename)[0]}.trimlog']
        
        # The funcion returns the names of the trimmed fastq.gz files
        return_value = [f'{outdir}/trim_{file_basename}']
        # return_value = [value[2:] if value[:2] == './' else value for value in return_value]
    
    # Call bbduk.sh
    subprocess.run(['bbduk.sh'] + arguments)

    # Return output trimmed fastq.gz file names
    return return_value
