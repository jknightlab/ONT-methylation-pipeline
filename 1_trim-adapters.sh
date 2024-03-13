#!/usr/bin/bash

## Author:      Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission

#SBATCH -A jknight.prj 
#SBATCH -J porechop

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/porechop.%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/porechop.%j.err 

#SBATCH -p short 
#SBATCH -c 1
##########################################################################################

# Set default parameter values
input_dir=$PWD
output_dir=$PWD
porechop="/well/jknight/users/awo868/software/Porechop/porechop-runner.py"
samtools="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/samtools/bin/samtools"

# Read in arguments
while getopts i:r:w:s:o:g:h opt
do
	case $opt in
	s)
		sample_list_file=$OPTARG
		;;
	i)
		input_dir=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	h)
		echo "Usage:	trim-reads.sh [-s sample_list] [-i input_dir] [-o output_dir]"
		echo ""
		echo "Where:"
		echo "-s		Text file containing a list of sample names (one per line) to be aligned. These names should match the naming convention of fastq files"
		echo "-i		Directory where input fastq files are located [defaults to the working directory]"
		echo "-o		Directory where output trimmed fastq files will be written [defaults to the working directory]"
		echo ""
		exit 1
		;;
	esac
done

# Output relevant information on how the job was run
echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} "
echo "------------------------------------------------" 

echo "[trim-reads.sh]:	Validating arguments..."

if [[ ! -f $sample_list_file ]]
then
        echo "[trim-reads.sh]:	ERROR: Sample list file not found."
        exit 2
fi 

if [[ ! -d $input_dir ]]
then
        echo "[trim-reads.sh]:	ERROR: Input (fastq) directory not found."
        exit 2
fi 

if [[ ! -d $output_dir ]]
then
        echo "[trim-reads.sh]:	ERROR: Output directory not found."
        exit 2
fi 

echo "[trim-reads.sh]:	Loading required modules..."
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"

# Finding sample of interest in sample list
echo "[trim-reads.sh]:	Reading in sample list..."
readarray sample_list < $sample_list_file
sample_name=$(echo ${sample_list[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

# Creating temporary output directory
if [[ ! -d "${output_dir}/tmp" ]]
then
        mkdir "${output_dir}/tmp"
fi

# Converting reads to FASTQ format
$samtools fastq -T "*" ${input_dir}/${sample_name}.bam > ${output_dir}/tmp/${sample_name}.fastq

# Trimming reads with Porechop
echo "[trim-reads.sh]:	Trimming adapter sequences for sample ${sample_name}..."
$porechop \
	-i "${output_dir}/tmp/${sample_name}.fastq" \
	-o "${output_dir}/${sample_name}_trimmed.fastq"

echo "[trin-reads.sh]: Compressing files..."
gzip "${output_dir}/${sample_name}_trimmed.fastq"

echo "[trim-reads.sh]:	Cleaning up..."
rm ${output_dir}/tmp/${sample_name}.fastq

echo "[trim-reads.sh]:	...done!"
