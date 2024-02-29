#!/usr/bin/bash

 ## Author:	Kiki Cano-Gamez (kiki.canogamez@well.ox.ac.uk)

##########################################################################################
# Specifying Slurm parameters for job submission
#SBATCH -A jknight.prj 
#SBATCH -J methylation-calls

#SBATCH -o /well/jknight/users/awo868/logs/ONT-pipeline/get-modified-base-calls_%j.out 
#SBATCH -e /well/jknight/users/awo868/logs/ONT-pipeline/get-modified-base-calls_%j.err 

#SBATCH -p long 
#SBATCH -c 2
##########################################################################################

# Setting default parameter values
input_dir=$PWD
output_dir=$PWD
modification_type="5mC"
reference_genome="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/resources/genome-references/minimap2/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
base_context="cpg"

# Specifying software paths
samtools="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/samtools/bin/samtools"
modbam2bed="/well/jknight/projects/sepsis-immunomics/cfDNA-methylation/ONT/software/modbam2bed/modbam2bed"

# Reading in arguments
while getopts i:o:s:m:c:r:h opt
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
	m)
		modification_type=$OPTARG
		;;
	c)
		base_context=$OPTARG
		;;
	r)
		reference_genome=$OPTARG
		;;
	h)
		echo "Usage:	get-modified-base-calls.sh [-i input_dir] [-o output_dir] [-s sample_list] [-c base_context] [-m modification_type] [-r reference_genome]"
		echo ""
		echo "Where:"
		echo "-i		Path to input directory containing aligned ONT reads (in BAM format). [defaults to the working directory]"
		echo "-o		Path to output directory where to write output BED files with base modification proportions per site [defaults to the working directory]"
		echo "-s		Path to a text file containing a list of samples (one sample per line). Sample names should match file naming patterns."
		echo "-c		Sequence context where to test for the presence of modified bases (either chg, chh or cpg). [defaults to cpg]"
		echo "-m		Type of modification to analyse (e.g. 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao, modA, moxC, modG, modT, modU and modN. For more information see modbam2bed's documentation) [defaults to 5mC]"
		echo "-r		Path to a reference genome FASTA file used for alignment. [defaults to a local GRCh38 reference genome FASTA file]"
		echo ""
		exit 1
		;;
	esac
done


# Validating arguments
echo "[get-modification-calls]:	Validating arguments..."

if [[ ! -d $input_dir ]]
then
		echo "[get-modification-calls]:	ERROR: Input directory not found."
		exit 2
fi 

if [[ ! -d $output_dir ]]
then
		echo "[get-modification-calls]:	ERROR: Output directory not found."
        exit 2
fi 

if [[ ! -f $sample_list ]]
then
        echo "[get-modification-calls]:	ERROR: Sample list file not found"
        exit 2
fi

if [[ ! -f $reference_genome ]]
then
        echo "[get-modification-calls]:	ERROR: Reference genome (FASTA) file not found"
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
echo "[get-modification-calls]:	Reading sample list..."
readarray sampleList < $sample_list

# Per-task processing 
## Defining input sample and barcode names
sampleName=$(echo ${sampleList[$((${SLURM_ARRAY_TASK_ID}-1))]} | sed 's/\n//g')

echo "[get-modification-calls]:	Setting up output directory structure..."
if [[ ! -d "${input_dir}/tmp" ]]
then
        mkdir "${input_dir}/tmp"
fi

if [[ ! -d "${output_dir}/${modification_type}" ]]
then
        mkdir "${output_dir}/${modification_type}"
fi

if [[ ! -d "${output_dir}/${modification_type}/${sampleName}" ]]
then
        mkdir "${output_dir}/${modification_type}/${sampleName}"
fi
 
## Indexing BAM files
echo "[get-modification-calls]:	Sorting BAM file ($sampleName)..."
$samtools sort  ${input_dir}/${sampleName}_aligned-reads.bam -o ${input_dir}/tmp/${sampleName}_aligned-reads_sorted.bam

echo "[get-modification-calls]:	Indexing BAM file ($sampleName)..."
$samtools index ${input_dir}/tmp/${sampleName}_aligned-reads_sorted.bam

echo "[get-modification-calls]:	Fetching modified base information and converting to BED format ($sampleName)..."
cd ${output_dir}/${modification_type}/${sampleName}
$modbam2bed \
	-m $modification_type \
	--aggregate \
	--$base_context \
	$reference_genome \
	${input_dir}/tmp/${sampleName}_aligned-reads_sorted.bam > \
	${sampleName}_modified-base-calls_${modification_type}_per-strand.bed

echo "[get-modification-calls]:	...done!"

