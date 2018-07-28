#!/bin/bash

#SBATCH -p Lewis
#SBATCH --account=biocommunity
#SBATCH -J GVCFgeno
#SBATCH --mem 200G
#SBATCH -N1
#SBATCH -n20
#SBATCH -t 2-00:00
#SBATCH --output=gtGVCFgeno-%A_%a-%j.out

## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send

#------------------------------------------------------------------
#USER defined option

# Path to ref genome
REFPATH="/home/buckleyrm/storage.lyonslab/cat_ref/"
REFNAME="Felis_catus_9.0.fa"
# Path to where gvcfs are kept
GVCFPATH="/home/buckleyrm/storage.lyonslab/results/Felis_catus/gvcf/"
# file containg a list of lab ids
LISTPATH="/home/buckleyrm/storage.lyonslab/dom_cat_run/"
LISTNAME="id.list"
# name of the output vcf
OUTPATH="/storage/hpc/group/UMAG/WORKING/buckleyrm/gvcf_geno/"
OUTNAME="dom_cat_run_186"
#------------------------------------------------------------------

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

# invoke with sbatch --array=1-$(ls ~/storage.lyonslab/cat_ref/target_loci/ | wc -l)%2 GVCFgeno.sh 

TARGETS=$(ls $REFPATH/target_loci/)
TARGET=$(echo $TARGETS | cut -d " " -f $SLURM_ARRAY_TASK_ID)


LEN=$(wc -l $LISTPATH/$LISTNAME | cut -f1 -d" ")
END=$(seq $LEN $(expr $LEN / -19) 1)
START=$(echo $(seq $(expr $LEN - $(expr $LEN / 19) + 1) $(expr $LEN / -19) 1) 1)


# clean dir

if [ -f $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list ]; then
	rm $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list
fi


for i in $(seq 1 20); do

	(
	sleep $((RANDOM % 20))
	
	if [ -f $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohort_$i.g.vcf.gz ]; then
                rm $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohort_$i.g.vcf.gz
        fi

	awk -v start=$(echo $START | cut -f $i -d " ") -v end=$(echo $END | cut -f $i -d " ") 'NR >= start && NR <= end { print }' $LISTPATH/$LISTNAME | 
	awk -v gvcf=$GVCFPATH '{print gvcf "/" $0 ".g.vcf.gz"}' > $OUTPATH/tmp.${TARGET%\.intervals}.cohort_$i.$LISTNAME

	echo $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohort_$i.g.vcf.gz &>> $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list

	java -Djava.io.tmpdir=$GVCFPATH/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
	-T CombineGVCFs \
	-R $REFPATH/$REFNAME \
	-V $OUTPATH/tmp.${TARGET%\.intervals}.cohort_$i.$LISTNAME \
	-L $REFPATH/target_loci/$TARGET \
	--log_to_file $(pwd)/gtCombine.${TARGET%\.intervals}.cohort_$i.log \
	--out $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohort_$i.g.vcf.gz

	#rm $OUTPATH/tmp.${TARGET%\.intervals}.cohort_$i.$LISTNAME
	)&

done

wait


java -Djava.io.tmpdir=$GVCFPATH/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nt 20 \
-T GenotypeGVCFs \
-R $REFPATH/$REFNAME \
-V $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list \
-L $REFPATH/target_loci/$TARGET \
--log_to_file $(pwd)/gtGenotype.${TARGET%\.intervals}.log
--out $OUTPATH/$OUTNAME.${TARGET%\.intervals}.vcf.gz

#cat $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list | xargs rm
#rm $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list



