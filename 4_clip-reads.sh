#!/usr/bin/bash

## Author:      Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission

#SBATCH -A jknight.prj 
#SBATCH -J clip-reads

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/clip-reads.%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/clip-reads.%j.err 

#SBATCH -p short 
#SBATCH -c 2
##########################################################################################

# Set default parameter values
input_dir=$PWD
output_dir=$PWD

modkit="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/modkit/modkit"

# Read in arguments
while getopts s:i:o:h opt
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
		echo "Usage:	clip-reads.sh [-s sample_list] [-i input_dir] [-o output_dir]"
		echo ""
		echo "Where:"
		echo "-s		Text file containing a list of sample names (one per line) to be aligned. These names should match the naming convention of BAM files"
		echo "-i		Directory where input trimmed, aligned, and MM tag-repaired reads (in BAM format) are located [defaults to the working directory]"
		echo "-o		Directory where output clipped BAM files will be written [defaults to the working directory]"
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

echo "[clip-reads.sh]:	Validating arguments..."

if [[ ! -f $sample_list_file ]]
then
        echo "[clip-reads.sh]:	ERROR: Sample list file not found."
        exit 2
fi 

if [[ ! -d $input_dir ]]
then
        echo "[clip-reads.sh]:	ERROR: Input directory not found."
        exit 2
fi 

if [[ ! -d $output_dir ]]
then
        echo "[clip-reads.sh]:	ERROR: Output directory not found."
        exit 2
fi 

# Parallelising task by sample
echo "[clipe-reads.sh]:	Reading in sample list..."
readarray sample_list < $sample_list_file
sample_name=$(echo ${sample_list[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

# Clipping methylation information at read ends
echo "[clip-reads.sh]:	Clipping methylation at read ends with modkit (${sample_name})..."
$modkit adjust-mods \
	--edge-filter 5,30 \
	${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired.bam \
	${output_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped.bam
	
echo "[clip-reads.sh]:	...done!"

