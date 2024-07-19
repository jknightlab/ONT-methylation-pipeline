#!/usr/bin/bash

## Author:      Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission

#SBATCH -A jknight.prj 
#SBATCH -J nanomix

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/nanomix.%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/nanomix.%j.err 

#SBATCH -p short 
#SBATCH -c 1
##########################################################################################

# Set default parameter values
input_dir=$PWD
output_dir=$PWD

methylation_atlas='/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/nanomix/atlases/39Bisulfite.tsv'
methylation_context='5mC'

# Read in arguments
while getopts s:i:o:a:c:h opt
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
	a)
		methylation_atlas=$OPTARG
		;;
	c)
		methylation_context=$OPTARG
		;;
	h)
		echo "Usage:	deconvolute-with-nanomix.sh [-s sample_list] [-i input_dir] [-o output_dir] [-a methylation_atlas] [-m methylation_context]"
		echo ""
		echo "Where:"
		echo "-s		Text file containing a list of sample names (one per line) to be deconvoluted. These names should match the naming of methylome files"
		echo "-i		Directory containing methylation pile up files for each sample (in BED format, as outputed by modkit). [defaults to the working directory]"
		echo "-o		Directory where output files will be written [defaults to the working directory]"
		echo "-a		Text file (TSV format) containing a methylation reference atlas to be used for deconvolution [defaults to a local file with the cell type atlas publised by Loyfer et al.]"
		echo "-c		Context in which to look at methylation (either '5mC' or '5hmC') [defaults to '5mC']"
		echo ""
		exit 1
		;;
	esac
done


# Output relevant information on how the job will be run
echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "Executing task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID} "
echo "------------------------------------------------" 

echo "[nanomix.sh]:	Validating arguments..."

if [[ ! -f $sample_list_file ]]
then
        echo "[nanomix.sh]:	ERROR: Sample list file not found."
        exit 2
fi 

if [[ ! -f $methylation_atlas ]]
then
        echo "[nanomix.sh]:	ERROR: Methylation reference atlas not found."
        exit 2
fi 

if [[ ! -d $input_dir ]]
then
        echo "[nanomix.sh]:	ERROR: Input directory not found."
        exit 2
fi 

if [[ ! -d $output_dir ]]
then
        echo "[nanomix.sh]:	ERROR: Output directory not found."
        exit 2
fi 

if [[ $methylation_context != "5mC" && $methylation_context != "5hmC" ]]
then
        echo "[nanomix.sh]:	ERROR: Methylation context is invalid. Please use either '5mC' or '5hmC'."
        exit 2
fi 

# Loading required modules and virtual environments
echo "[nanomix.sh]:	Loading required software modules..."
module load Python/3.10.8-GCCcore-12.2.0
module load Rust/1.65.0-GCCcore-12.2.0

echo "[nanomix.sh]:	Activating virtual environment..."
source /well/jknight/users/awo868/python/nanomix-${MODULE_CPU_TYPE}/bin/activate

# Parallelising task by sample
echo "[nanomix.sh]:	Reading in sample list..."
readarray sample_list < $sample_list_file
sample_name=$(echo ${sample_list[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

# Creating methylome files
echo "[nanomix.sh]:	Creating methylome files in nanomix format (${sample_name})..."
if [[ ! -d "${output_dir}/tmp" ]]
then
		mkdir "${output_dir}/tmp"
fi

if [[ $methylation_context == "5mC" ]] 
then
        cat ${input_dir}/${sample_name}_methylation-pileup.bed | \
        	awk '$4=="m"' | \
        	awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$12}' > \
        	${output_dir}/tmp/${sample_name}_5mC_methylome.bed
        
        echo -e "chr\tstart\tend\ttotal_calls\tmodified_calls" | \
			cat - ${output_dir}/tmp/${sample_name}_5mC_methylome.bed > \
			${output_dir}/tmp/${sample_name}_5mC_methylome.tsv
fi

if [[ $methylation_context == "5hmC" ]] 
then
		cat ${input_dir}/${sample_name}_methylation-pileup.bed | \
			awk '$4=="h"' | \
			awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$12}' > \
			${output_dir}/tmp/${sample_name}_5hmC_methylome_no-header.bed
			
		echo -e "chr\tstart\tend\ttotal_calls\tmodified_calls" | \
			cat - ${output_dir}/tmp/${sample_name}_5hmC_methylome.bed > \
			${output_dir}/tmp/${sample_name}_5hmC_methylome.tsv
fi


# Deconvoluting methylome
echo "[nanomix.sh]:	Deconvoluting methylome with nanomix (${sample_name})..."
if [[ $methylation_context == "5mC" ]] 
then
		nanomix deconvolute \
			-a $methylation_atlas \
			${output_dir}/tmp/${sample_name}_5mC_methylome.tsv > \
			${output_dir}/${sample_name}_tissue-proportions_nanomix_5mC.txt
fi

if [[ $methylation_context == "5hmC" ]] 
then
		nanomix deconvolute \
			-a $methylation_atlas \
			${output_dir}/tmp/${sample_name}_5hmC_methylome.tsv > \
			${output_dir}/${sample_name}_tissue-proportions_nanomix_5hmC.txt
fi

# Removing temporary files
echo "[nanomix.sh]:	Cleaning up..."	
if [[ $methylation_context == "5mC" ]] 
then
		rm "${output_dir}/tmp/${sample_name}_5mC_methylome.bed"
		rm "${output_dir}/tmp/${sample_name}_5mC_methylome.tsv"
fi

if [[ $methylation_context == "5hmC" ]] 
then
		rm "${output_dir}/tmp/${sample_name}_5hmC_methylome.bed"
		rm "${output_dir}/tmp/${sample_name}_5hmC_methylome.tsv"
fi

echo "[nanomix.sh]: ...done!"

