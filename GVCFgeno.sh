#!/bin/bash

#SBATCH -p BioCompute
#SBATCH -J GVCFgeno
#SBATCH --mem 50G
#SBATCH -N1
#SBATCH -n28
#SBATCH -t 0-12:00
#SBATCH -o gtGVCFgeno.o%j
#SBATCH -e gtGVCFgeno.e%j

## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send

#------------------------------------------------------------------
#USER defined option

# Path to ref genome
REFPATH="/home/buckleyrm/storage.lyonslab/cat_ref/"
REFNAME="Felis_catus_9.0.fa"
# Path to where gvcfs are kept
GVCFPATH="/home/buckleyrm/storage.lyonslab/big.test.preprocess/gvcf"
# file containg a list of lab ids
LISTPATH="/home/buckleyrm/storage.lyonslab/big.test.preprocess/vcf/"
LISTNAME="id.list"
# name of the output vcf 
OUTPATH="/home/buckleyrm/storage.lyonslab/big.test.preprocess/vcf"
OUTNAME="big.test"
#------------------------------------------------------------------

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

awk -v gvcf=$GVCFPATH '{print gvcf "/" $0 ".sort.markDup.realign.bam.g.vcf.gz"}' $LISTPATH/$LISTNAME > $LISTPATH/tmp.$LISTNAME

java -Djava.io.tmpdir=$GVCFPATH/tmp -XX:ParallelGCThreads=2 -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nt 26 \
-T GenotypeGVCFs \
-R $REFPATH/$REFNAME \
-V $LISTPATH/tmp.$LISTNAME \
--out $OUTPATH/$OUTNAME.vcf.gz


rm $LISTPATH/tmp.$LISTNAME
