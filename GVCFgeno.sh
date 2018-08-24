#!/bin/bash

#SBATCH -p BioCompute
#SBATCH --account=biocommunity
#SBATCH -J GVCFgeno
#SBATCH --mem 250G
#SBATCH -N1
#SBATCH -n26
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
OUTNAME="dom_cat_run_183"
#------------------------------------------------------------------

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

# invoke with sbatch --array=1-$(ls ~/storage.lyonslab/cat_ref/target_loci/ | wc -l)%2 GVCFgeno.sh 


THREADS=20



TARGETS=$(ls $REFPATH/target_loci/)
TARGET=$(echo $TARGETS | cut -d " " -f $SLURM_ARRAY_TASK_ID)


LEN=$(wc -l $LISTPATH/$LISTNAME | cut -f1 -d" ")

if [ $(eval $LEN % $THREADS) -eq 0 ]; then
	END=$(seq $LEN -$(expr $LEN / $THREADS) 1)
	START=$(seq $(expr $LEN - $(expr $LEN / $THREADS) + 1) -$(expr $LEN / $THREADS) 1)
elif [ $(eval $LEN % $THREADS) -eq 1 ]; then
	END=$(seq $LEN -$(expr $LEN / $THREADS) $(expr $LEN / $THREADS))
	START=$(echo $(seq $(expr $LEN - $(expr $LEN / $THREADS) + 1) -$(expr $LEN / $THREADS) $(expr $LEN / $THREADS)) 1)
else
	END=$(seq $LEN -$(expr $(expr $LEN / $THREADS) + 1) 1)
	START=$(echo $(seq $(expr $LEN - $(expr $LEN / $THREADS)) -$(expr $(expr $LEN / $THREADS) + 1) 1) 1)
fi


# clean dir

if [ -f $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list ]; then
	rm $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list
fi


for i in $(seq 1 $(echo $END | wc -w)); do

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

	rm $OUTPATH/tmp.${TARGET%\.intervals}.cohort_$i.$LISTNAME
	)&

done

wait

for i in $(seq 1 $(echo $END | wc -w)); do
	if [ -f $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohort_$i.g.vcf.gz.tbi ];
	then
		echo ${TARGET%\.intervals} cohort_$i is done, continuing
	else
		echo ${TARGET%\.intervals} cohort_$i is not done, missing .tbi file
		wait
	fi
done

java -Djava.io.tmpdir=$GVCFPATH/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nt 26 \
-T GenotypeGVCFs \
-R $REFPATH/$REFNAME \
-V $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list \
-L $REFPATH/target_loci/$TARGET \
--log_to_file $(pwd)/gtGenotype.${TARGET%\.intervals}.log \
--out $OUTPATH/$OUTNAME.${TARGET%\.intervals}.vcf.gz

cat $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list | xargs rm
rm $OUTPATH/$OUTNAME.${TARGET%\.intervals}.cohorts.list



