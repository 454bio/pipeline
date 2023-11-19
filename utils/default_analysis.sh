#!/bin/bash

# default paths for binaries and tools
basecall_path=~/Code/pipeline/basecall
util_path=~/Code/pipeline/utils
align_path=~/Code/pipeline/align

# defaults for reference genome
ref_path=~/Code/pipeline/refs
ref_name=synth_GC23.fasta
#ref_name=PhiX_GC23.fasta

# we require an input run name, the output files generated will all begin with this name
run_name=$1
if [ -z "$1" ]
then
    echo missing run name, using default 'run'
    run_name=run
fi

options=none
if [[ "$2" == "fast" ]]
then
    options=fast
fi


if [[ "${options}" == "fast" ]]
then
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.15 0.11 0.017 > ${run_name}.log
else
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -g > ${run_name}.log
fi

python3 ${util_path}/fastq_stats.py -i ${run_name}.fastq >> ${run_name}.log
python3 ${align_path}/align_reads.py --ref ${ref_path}/${ref_name} --in ${run_name}.fastq --out ${run_name}.txt --poslist 0,16,32,48 --direction 1 >> ${run_name}.log
#python3 ${align_path}/align_reads.py --ref ${ref_path}/${ref_name} --in ${run_name}.fastq --out ${run_name}.txt >> ${run_name}.log
python3 ${util_path}/filtered_stats.py -i ${run_name}.txt.filtered >> ${run_name}.log
python3 ${util_path}/report.py -i ${run_name}.log -o ${run_name}.html

