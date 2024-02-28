#!/usr/bin/bash

## Author:	Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission
#SBATCH -A jknight.prj 
#SBATCH -J trim-adapteres

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/trim-adapters_%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/trim-adapters_%j.err 

#SBATCH -p short 
#SBATCH -c 1
##########################################################################################

# Setting default parameter values
input_dir=$PWD
output_dir=$PWD


# Reading in arguments
while getopts i:o:s:b:h opt
do
	case $opt in
	i)
		input_dir=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	s)
		sample_list=$OPTARG
		;;
	h)
		echo "Usage:	trim-adapters.sh [-i input_dir] [-o output_dir] [-s sample_list]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing merged base calls obtained from Nanopore's Dorado/MinKnow software. Calls must be in BAM format and there should be one fil per sample. [defaults to the working directory]"
		echo "-o		Path to output directory where to write trimmed FASTQ files with base calls (these files will retained information tags regarding modified bases) [defaults to the working directory]"
		echo "-s		Path to a text file containing a list of samples (one sample per line). Sample names should match file naming patterns."
		echo ""
		exit 1
		;;
	esac
done


# Validating arguments
echo "[trim-adapters]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
		echo "[trim-adapters]:	ERROR: Input directory not found."
		exit 2
fi 

if [[ ! -d $output_dir ]]
then
		echo "[trim-adapters]:	ERROR: Output directory not found."
        exit 2
fi 

if [[ ! -f $sample_list ]]
then
        echo "[trim-adapters]:	ERROR: Sample list file not found"
        exit 2
fi


# Outputing relevant information on how the job was run
echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} "
echo "------------------------------------------------" 


# Loading required modules and files
echo "[trim-adapters]:	Loading required modules..."
module load dorado/0.5.1-foss-2022a-CUDA-11.7.0

echo "[trim-adapters]:	Reading sample list..."
readarray sampleList < $sample_list


# Per-task processing 
## Defining input sample and barcode names
sampleName=$(echo ${sampleList[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

echo "[trim-adapters]:	Trimming ONT adapters from reads for $sampleName..."
dorado trim ${input_dir}/${sampleName}.bam > ${output_dir}/${sampleName}_trimmed.bam

echo "[trim-adapters]:	...done!"

