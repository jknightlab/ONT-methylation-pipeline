#!/usr/bin/bash

## Author:	Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission
#SBATCH -A jknight.prj 
#SBATCH -J align-reads

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/align-reads_%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/align-reads_%j.err 

#SBATCH -p long 
#SBATCH -c 3
##########################################################################################

# Setting default parameter values
input_dir=$PWD
output_dir=$PWD
reference_genome="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/resources/genome-references/minimap2/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Specifying software paths
samtools="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/samtools/bin/samtools"
minimap2="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/minimap2/minimap2"

# Reading in arguments
while getopts i:o:s:r:h opt
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
	r)
		reference_genome=$OPTARG
		;;
	h)
		echo "Usage:	align-reads.sh [-i input_dir] [-o output_dir] [-s sample_list]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing adapter-trimmed ONT base calls (in FASTQ format). [defaults to the working directory]"
		echo "-o		Path to output directory where to write aligned BAM files [defaults to the working directory]"
		echo "-s		Path to a text file containing a list of samples (one sample per line). Sample names should match file naming patterns."
		echo "-r		Path to a reference genome FASTA file used for alignment. [defaults to a local GRCh38 reference genome file for minimap2]"
		echo ""
		exit 1
		;;
	esac
done


# Validating arguments
echo "[align-reads]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
		echo "[align-reads]:	ERROR: Input directory not found."
		exit 2
fi 

if [[ ! -d $output_dir ]]
then
		echo "[align-reads]:	ERROR: Output directory not found."
        exit 2
fi 

if [[ ! -f $sample_list ]]
then
        echo "[align-reads]:	ERROR: Sample list file not found"
        exit 2
fi

if [[ ! -f $reference_genome ]]
then
        echo "[align-reads]:	ERROR: Reference genome (FASTA) file not found"
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


# Loading input files
echo "[align-reads]:	Reading sample list..."
readarray sampleList < $sample_list


# Per-task processing 
## Defining input sample and barcode names
sample_name=$(echo ${sampleList[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

echo "[align-reads]:	Aligning ONT reads with minimap2 ($sample_name)..."
$minimap2 -ax map-ont $reference_genome ${input_dir}/${sample_name}_trimmed.fastq.gz -y | \
	$samtools view -q 10 -b -o ${output_dir}/${sample_name}_trimmed_aligned-reads.bam

echo "[align-reads]:	...done!"

