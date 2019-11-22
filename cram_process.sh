#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition Lewis,hpc5,hpc4,BioCompute  # for jobs < 2hrs try 'General'
#SBATCH -N1
#SBATCH -n40 # cores 
#SBATCH --mem 100G  # memory 
#SBATCH -t 2-00:00  # days-hours:minutes
#SBATCH --account=lyonslab  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=cram_file
#SBATCH --output=cram-%j.out  # %j is the unique jobID
#
## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#
## array options
#
#-------------------------------------------------------------------------------
module load samtools/samtools-1.9

CRAMDIR=/home/buckleyrm/storage.lyonslab/sra_download/ena_dl

for i in $(ls $CRAMDIR/*bam); do
	samtools split --threads 40 $i  
done

for i in $(ls *bam); do
	samtools sort -n --threads 40 -m 2500M $i > sort.$i
done


