#!/bin/bash

#SBATCH -p Lewis,BioCompute,hpc4,hpc5
#SBATCH --account=biocommunity
#SBATCH -J rename
#SBATCH --mem 50G
#SBATCH -N1
#SBATCH -n1
#SBATCH -t 2-00:00
#SBATCH -o merge_chr.o%j
#SBATCH -e merge_chr.e%j

## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send



module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
module load bcftools/bcftools-1.8


TMP=$(pwd)/tmp
NAMEFILE=/home/buckleyrm/storage.lyonslab/users/buckleyrm/name_lists/old_new_names.tsv
VCFPATH=/home/buckleyrm/storage.hpc.buckelyrm/gvcf_geno/domestics_190823/filtered
VCFIN=dom_filter_190823
VCFSFX=filtered.sort.vcf.gz
VCFOUT=domestics_190823.filtered.named
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa


if [ -s $VCFPATH/$VCFIN.$VCFSFX.tbi ]; then
        echo -e "$(date)\nconcat vcf index found, skipping CatVariants ...\n\n"
else
        echo -e "$(date)\nconcat vcf index NOT found, running CatVaraiants ...\n\n"

java -Djava.io.tmpdir=$TMP -cp /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
-R $REF \
-out $VCFPATH/$VCFIN.$VCFSFX \
-V $VCFPATH/$VCFIN.chrA1.$VCFSFX \
-V $VCFPATH/$VCFIN.chrA2.$VCFSFX \
-V $VCFPATH/$VCFIN.chrA3.$VCFSFX \
-V $VCFPATH/$VCFIN.chrB1.$VCFSFX \
-V $VCFPATH/$VCFIN.chrB2.$VCFSFX \
-V $VCFPATH/$VCFIN.chrB3.$VCFSFX \
-V $VCFPATH/$VCFIN.chrB4.$VCFSFX \
-V $VCFPATH/$VCFIN.chrC1.$VCFSFX \
-V $VCFPATH/$VCFIN.chrC2.$VCFSFX \
-V $VCFPATH/$VCFIN.chrD1.$VCFSFX \
-V $VCFPATH/$VCFIN.chrD2.$VCFSFX \
-V $VCFPATH/$VCFIN.chrD3.$VCFSFX \
-V $VCFPATH/$VCFIN.chrD4.$VCFSFX \
-V $VCFPATH/$VCFIN.chrE1.$VCFSFX \
-V $VCFPATH/$VCFIN.chrE2.$VCFSFX \
-V $VCFPATH/$VCFIN.chrE3.$VCFSFX \
-V $VCFPATH/$VCFIN.chrF1.$VCFSFX \
-V $VCFPATH/$VCFIN.chrF2.$VCFSFX \
-V $VCFPATH/$VCFIN.chrM.$VCFSFX \
-V $VCFPATH/$VCFIN.chrX.$VCFSFX \
-V $VCFPATH/$VCFIN.UNMAPPED.$VCFSFX \
-assumeSorted

fi

if [ -s $VCFPATH/$VCFOUT.vcf.gz.tbi ]; then
        echo -e "$(date)\nnamed vcf index found, skipping reheader ...\n\n"
else
        echo -e "$(date)\nnamed vcf index NOT found, running reheader ...\n\n"
	bcftools reheader -s $NAMEFILE $VCFPATH/$VCFIN.$VCFSFX > $VCFPATH/$VCFOUT.vcf.gz
	tabix $VCFPATH/$VCFOUT.vcf.gz

fi


