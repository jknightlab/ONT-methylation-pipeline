#!/usr/bin/bash

## Author:      Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission

#SBATCH -A jknight.prj 
#SBATCH -J tag-repair

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/repair-tags.%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/repair-tags.%j.err 

#SBATCH -p short 
#SBATCH -c 3
##########################################################################################

# Set default parameter values
input_dir_raw=$PWD
input_dir_aligned=$PWD
output_dir=$PWD

samtools="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/samtools/bin/samtools"
modkit="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/modkit/modkit"

# Read in arguments
while getopts s:i:I:o:h opt
do
	case $opt in
	s)
		sample_list_file=$OPTARG
		;;
	i)
		input_dir_aligned=$OPTARG
		;;
	I)
		input_dir_raw=$OPTARG
		;;
	o)
		output_dir=$OPTARG
		;;
	h)
		echo "Usage:	repair-MM-tags.sh [-s sample_list] [-f input_dir_fastq] [-b input_dir_bam] [-o output_dir]"
		echo ""
		echo "Where:"
		echo "-s		Text file containing a list of sample names (one per line) to be aligned. These names should match the naming convention of fastq files"
		echo "-i		Directory where input trimmed and aligned reads (in BAM format) are located [defaults to the working directory]"
		echo "-I		Directory where the original untrimmed BAM files are located (these files must contain intact MM tags) [defaults to the working directory]"
		echo "-o		Directory where output tag-repaired BAM files will be written [defaults to the working directory]"
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

echo "[repair-tags.sh]:	Validating arguments..."

if [[ ! -f $sample_list_file ]]
then
        echo "[repair-tags.sh]:	ERROR: Sample list file not found."
        exit 2
fi 

if [[ ! -d $input_dir_aligned ]]
then
        echo "[repair-tags.sh]:	ERROR: Input (aligned BAM) directory not found."
        exit 2
fi 

if [[ ! -d $input_dir_raw ]]
then
        echo "[repair-tags.sh]:	ERROR: Input (raw BAM) directory not found."
        exit 2
fi 

if [[ ! -d $output_dir ]]
then
        echo "[repair-tags.sh]:	ERROR: Output directory not found."
        exit 2
fi 

# Parallelising task by sample
echo "[repair-tags.sh]:	Reading in sample list..."
readarray sample_list < $sample_list_file
sample_name=$(echo ${sample_list[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

# Creating temporary output directory
if [[ ! -d "${output_dir}/tmp" ]]
then
        mkdir "${output_dir}/tmp"
fi

# Sorting BAM files
echo "[repair-tags.sh]:	Sorting BAM files by read name ($sample_name)..."
$samtools sort -n -O bam -o ${output_dir}/tmp/${sample_name}_tag-acceptor.bam ${input_dir_aligned}/${sample_name}_trimmed_aligned-reads.bam
$samtools sort -n -O bam -o ${output_dir}/tmp/${sample_name}_tag-donor.bam ${input_dir_raw}/${sample_name}.bam

# Repairing MM tags
echo "[repair-tags.sh]:	Repairing MM tags with modkit (${sample_name})..."
$modkit repair \
	--donor-bam ${output_dir}/tmp/${sample_name}_tag-donor.bam \
	--acceptor-bam ${output_dir}/tmp/${sample_name}_tag-acceptor.bam \
	--log-filepath ${output_dir}/${sample_name}_modkit-repair.log \
	--output-bam ${output_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired.bam
	
echo "[repair-tags.sh]:	Cleaning up..."
rm "${output_dir}/tmp/${sample_name}_tag-acceptor.bam"
rm "${output_dir}/tmp/${sample_name}_tag-donor.bam"

echo "[repair-tags.sh]:	...done!"

