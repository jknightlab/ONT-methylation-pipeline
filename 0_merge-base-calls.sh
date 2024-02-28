#!/usr/bin/bash

## Author:	Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission
#SBATCH -A jknight.prj 
#SBATCH -J merge-base-calls

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/merge-calls_%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/merge-calls_%j.err 

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
	b)
		barcode_list=$OPTARG
		;;
	h)
		echo "Usage:	merge-calls.sh [-i input_dir] [-o output_dir] [-s sample_list] [-b barcode_list]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing base calls obtained from Nanopore's Dorado/MinKnow software. This directory should contain one subdirectory per barcode, with calls in BAM format contained within it. This is must have the same structure as the bam_pass directory outputed by the MinKnow software [defaults to the working directory]"
		echo "-o		Path to output directory where to write merged BAM files with base calls  [defaults to the working directory]"
		echo "-s		Path to a text file containing a list of samples (one sample per line). Sample names should match file naming patterns."
		echo "-b		Path to a text file containing a list of the Nanopore native barcodes used for each sample. Barcode order should match the corresponding order of the sample list file."
		echo ""
		exit 1
		;;
	esac
done

# Validating arguments
echo "[merge-calls]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
		echo "[merge-calls]:	ERROR: Input directory not found."
		exit 2
fi 

if [[ ! -d $output_dir ]]
then
		echo "[merge-calls]:	ERROR: Output directory not found."
        exit 2
fi 

if [[ ! -f $sample_list ]]
then
        echo "[merge-calls]:	ERROR: Sample list file not found"
        exit 2
fi


if [[ ! -f $barcode_list ]]
then
        echo "[merge-calls]:	ERROR: Barcode list file not found"
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
echo "[merge-calls]:	Loading required modules..."
module load samtools/1.8-gcc5.4.0

echo "[merge-calls]:	Reading sample list..."
readarray sampleList < $sample_list

echo "[merge-calls]:	Reading barcode list..."
readarray barcodeList < $barcode_list

# Per-task processing 
## Defining input sample and barcode names
sampleName=$(echo ${sampleList[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')
barcodeID=$(echo ${barcodeList[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

echo "[merge-calls]:	Merging base calls for sample $sampleName..."
samtools merge \
	${output_dir}/${sampleName}.bam \
	${input_dir}/${barcodeID}/*.bam
	
echo "[merge-calls]:	...done!"
