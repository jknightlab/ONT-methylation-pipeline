#!/usr/bin/bash

## Author:      Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission

#SBATCH -A jknight.prj 
#SBATCH -J methylation-calling

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/methylation-calling.%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/methylation-calling.%j.err 

#SBATCH -p short 
#SBATCH -c 3
##########################################################################################

# Setting default parameter values
input_dir=$PWD
output_dir=$PWD
reference_genome="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/resources/genome-references/minimap2/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Specifying software paths
samtools="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/samtools/bin/samtools"
modkit="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/modkit/modkit"

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
		echo "Usage:	get-methylation-calls.sh [-i input_dir] [-o output_dir] [-s sample_list] [-c base_context] [-m modification_type] [-r reference_genome]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing aligned ONT reads (in BAM format). [defaults to the working directory]"
		echo "-o		Path to output directory where to write output BED files with base modification pileups and prpotions per site [defaults to the working directory]"
		echo "-s		Path to a text file containing a list of samples (one sample per line). Sample names should match file naming patterns."
		echo "-r		Path to a reference genome FASTA file used for alignment. [defaults to a local GRCh38 reference genome FASTA file]"
		echo ""
		exit 1
		;;
	esac
done

# Validating arguments
echo "[methylation-calls]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
		echo "[methylation-calls]:	ERROR: Input directory not found."
		exit 2
fi 

if [[ ! -d $output_dir ]]
then
		echo "[methylation-calls]:	ERROR: Output directory not found."
        exit 2
fi 

if [[ ! -f $sample_list ]]
then
        echo "[methylation-calls]:	ERROR: Sample list file not found"
        exit 2
fi

if [[ ! -f $reference_genome ]]
then
        echo "[methylation-calls]:	ERROR: Reference genome (FASTA) file not found"
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


# Parallelising task by sample
echo "[methylation-calls]:	Reading in sample list..."
readarray sampleList < $sample_list
sample_name=$(echo ${sampleList[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

# Indexing and sorting BAM files
if [[ ! -f "${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped_sorted.bam" ]]
then
	echo "[methylation-calls]:	Sorting BAM file ($sample_name)..."
	$samtools sort ${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped.bam -o ${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped_sorted.bam
fi

if [[ ! -f "${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped_sorted.bam.bai" ]]
then
	echo "[methylation-calls]:	Indexing BAM file ($sample_name)..."
	$samtools index ${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped_sorted.bam
fi

# Creating pileup BED files
echo "[methylation-calls]:	Creating methylation BED files with modkit ($sample_name)..."
$modkit pileup \
	${input_dir}/${sample_name}_trimmed_aligned-reads_tag-repaired_clipped_sorted.bam \
	${output_dir}/${sample_name}_methylation-pileup.bed \
	--cpg \
	--ref $reference_genome

echo "[methylation-calls]:	...done!"

