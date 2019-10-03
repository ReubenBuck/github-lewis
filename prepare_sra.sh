#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p BioCompute,Lewis,hpc5
#SBATCH --account=lyonslab
#SBATCH -J BQSR_Strict
#SBATCH --output=sra_bam_prep.out
#SBATCH --mem 50G
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 1-00:00


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send


module load samtools/samtools-1.9-test
module load pigz/pigz-2.4

#------------------------------------------------------
# USER CONFIG
#------------------------------------------------------

CWD=/home/buckleyrm/storage.lyonslab/users/buckleyrm/test_bam2fastq
TMP=felCat.Fcat22055.Chediak


#-----------------------------------------------------

#invoke with:
# BAMLIST=<bamPathList> sbatch prepare_sra.sh

if [ ! -d $CWD/$TMP ]; then 
	mkdir -p $CWD/$TMP
fi

for bam in $(cat $BAMLIST); do
	bamFile=$(echo $bam | sed "s|^.*/||")
	samtools fastq -@ 19 -1 $CWD/$TMP/${bamFile}.R1.fastq -2 $CWD/$TMP/${bamFile}.R2.fastq $bam
done

cat $CWD/$TMP/*.R1.fastq >> $CWD/$TMP.R1.fastq &
cat $CWD/$TMP/*.R2.fastq >> $CWD/$TMP.R2.fastq &

wait

pigz -p 20 $CWD/$TMP.R1.fastq
pigz -p 20 $CWD/$TMP.R2.fastq






