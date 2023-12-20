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

options=$2


if [[ "${options}" == "fast" ]]; then
    echo running fast
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.16 0.11 0.017 > ${run_name}.log
elif [[ "${options}" == "gridw" ]]; then
    echo running grid with custom params and whitelist
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -g -v -G 0.20 0.30 21 0.10 0.17 15 0.02 0.07 11 -w ${run_name}.whitelist.txt -c 7  > ${run_name}.log
elif [[ "${options}" == "fastopt" ]]; then
    echo running fast optimized with normalzed data
    #${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.27 0.145 0.05 > ${run_name}.log
    #${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.29 0.14 0.045 > ${run_name}.log
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.22 0.135 0.045 > ${run_name}.log
elif [[ "${options}" == "fastn" ]]; then
    echo running fast with normalzed data
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.15 0.11 0.017 -n > ${run_name}.log
elif [[ "${options}" == "gridnw" ]]; then
    echo running grid with custom params and whitelist
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -n -g -v -G 0.20 0.30 21 0.10 0.17 15 0.02 0.07 11 -w ${run_name}.whitelist.txt -c 7  > ${run_name}.log
elif [[ "${options}" == "fastnopt" ]]; then
    echo running fast optimized with normalzed data
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.27 0.145 0.05 -n > ${run_name}.log
elif [[ "${options}" == "fastnw" ]]; then
    echo running fast with normalzed data and white-list
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -p 0.15 0.11 0.017 -n -w ${run_name}.whitelist.txt > ${run_name}.log
elif [[ "${options}" == "gridp" ]]; then
    echo running grid with custom params
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -n -g -v -G 0.10 0.20 11 0.05 0.15 11 0.02 0.05 4  > ${run_name}.log
else
    echo running grid search
    ${basecall_path}/Basecall -i color_transformed_spots.csv -o ${run_name}.fastq -g > ${run_name}.log
fi

python3 ${util_path}/fastq_stats.py -i ${run_name}.fastq >> ${run_name}.log
python3 ${align_path}/align_reads.py --ref ${ref_path}/${ref_name} --in ${run_name}.fastq --out ${run_name}.txt --poslist 0,16,32,48 --direction 1 >> ${run_name}.log
#python3 ${align_path}/align_reads.py --ref ${ref_path}/${ref_name} --in ${run_name}.fastq --out ${run_name}.txt >> ${run_name}.log
python3 ${util_path}/filtered_stats.py -i ${run_name}.txt.filtered >> ${run_name}.log
python3 ${util_path}/report.py -i ${run_name}.log -o ${run_name}.html

