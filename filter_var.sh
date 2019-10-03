#!/bin/bash
#SBATCH -p BioCompute,htc4,hpc5,Lewis
#SBATCH --account=lyonslab
#SBATCH -J Filter
#SBATCH --mem 50G
#SBATCH -N1
#SBATCH -n5
#SBATCH -t 2-00:00
#SBATCH --output=Filter-%A_%a-%j.out

## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send

#------------------------------------------------------------------
#USER defined option
REFPATH=/home/buckleyrm/storage.lyonslab/cat_ref
REF=Felis_catus_9.0.fa

GATK=/cluster/software/gatk/gatk-3.8/
PICARD=/cluster/software/picard-tools/picard-tools-2.1.1

INDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/gvcf_geno/domestics_190823
PREFIX=domestics_190823

OUTDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/gvcf_geno/domestics_190823/filtered
OUT=dom_filter_190823

TMPDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/TMP
LOGDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/gvcf_geno/domestics_190823/filtered


#------------------------------------------------------------------

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
module load picard-tools/picard-tools-2.1.1

# invoke with sbatch --array=1-$(ls ~/storage.lyonslab/cat_ref/target_loci/ | wc -l) filter_var.sh 


TARGETS=$(ls $REFPATH/target_loci/ | sed "s/.intervals//g")
TARGET=$(echo $TARGETS | cut -d " " -f $SLURM_ARRAY_TASK_ID)


#extract snps
if [ ! -d $LOGDIR ]; then
	mkdir -p $LOGDIR
fi

if [ -f $LOGDIR/$OUT.$TARGET.run.log ]; then
	touch $LOGDIR/$OUT.$TARGET.run.log 
fi

if [ ! -d $TMPDIR/$OUT.$TARGET ]; then
        mkdir -p $TMPDIR/$OUT.$TARGET
fi


echo -e "Start on $(date)" &>> $LOGDIR/$OUT.$TARGET.run.log


if [ -s $OUTDIR/$OUT.$TARGET.snps.vcf.gz.tbi ]; then
	echo -e "$(date)\nraw SNP VCF index found, skipping selection ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
else
	echo -e "$(date)\nraw SNP VCF index NOT found, selecting ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
	(
		java -jar $GATK/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $REFPATH/$REF \
		-V $INDIR/$PREFIX.$TARGET.vcf.gz \
		-selectType SNP \
		--log_to_file $LOGDIR/$OUT.$TARGET.snps.log \
		-o $OUTDIR/$OUT.$TARGET.snps.vcf.gz 
		)&
fi



sleep 5s

if [ -s $OUTDIR/$OUT.$TARGET.indels.vcf.gz.tbi ]; then
	echo -e "$(date)\nraw INDEL VCF index found, skipping selection ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
else
	echo -e "$(date)\nraw SNP VCF index NOT found, selecting ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
	(
		java -jar $GATK/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $REFPATH/$REF \
		-V $INDIR/$PREFIX.$TARGET.vcf.gz \
		-selectType INDEL \
		--log_to_file $LOGDIR/$OUT.$TARGET.indels.log \
		-o $OUTDIR/$OUT.$TARGET.indels.vcf.gz 
		)&
fi


wait

# filter the snp call set
if [ -s $OUTDIR/$OUT.$TARGET.filtered.snps.vcf.gz.tbi ]; then
	echo -e "$(date)\nfiltered SNP VCF indexes found, skipping filtering ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
else
	echo -e "$(date)\nfiltered SNP VCF indexes NOT found, filtering ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
	(
		java -jar $GATK/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R $REFPATH/$REF \
		-V $OUTDIR/$OUT.$TARGET.snps.vcf.gz \
		--filterExpression "QD < 2.0" \
		--filterName "QD" \
		--filterExpression "FS > 60.0" \
		--filterName "FS" \
		--filterExpression "SOR > 3.0" \
		--filterName "SOR" \
		--filterExpression "ReadPosRankSum < -8.0" \
		--filterName "ReadPosRankSum" \
		--filterExpression "MQ < 40.0" \
		--filterName "MQ" \
		--filterExpression "MQRankSum < -12.5" \
		--filterName "MQRankSum" \
		--log_to_file $LOGDIR/$OUT.$TARGET.filtered.snps.log \
		-o $OUTDIR/$OUT.$TARGET.filtered.snps.vcf.gz 
		)&

fi

# filter the indels call set
if [ -s $OUTDIR/$OUT.$TARGET.filtered.indels.vcf.gz.tbi ]; then
	echo -e "$(date)\nfiltered INDEL VCF index found, skipping filtering ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
else

	echo -e "$(date)\nfiltered INDEL VCF indexes NOT found, filtering ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
	#GATK
	sleep 5s
	(
		java -jar $GATK/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R $REFPATH/$REF \
		-V $OUTDIR/$OUT.$TARGET.indels.vcf.gz \
		--filterExpression "QD < 2.0" \
		--filterName "QD" \
		--filterExpression "FS > 200.0" \
		--filterName "FS" \
		--filterExpression "SOR > 10.0" \
		--filterName "SOR" \
		--filterExpression "ReadPosRankSum < -20.0" \
		--filterName "ReadPosRankSum" \
		--log_to_file $LOGDIR/$OUT.$TARGET.filtered.indels.log \
		-o $OUTDIR/$OUT.$TARGET.filtered.indels.vcf.gz 
		)&

fi

wait


# create tmp dirs

if [[ ! -d $TMPDIR/$OUT.$TARGET ]]; then
	mkdir -p $TMPDIR/$OUT.$TARGET
fi

# sort vcf files

if [ -s $OUTDIR/$OUT.$TARGET.filtered.sort.vcf.gz.tbi ]; then
	echo -e "$(date)\nGATK filtered index files found, skipping sorting ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
else
	echo -e "$(date)\nGATK filtered index files NOT found, sorting ...\n\n" &>> $LOGDIR/$OUT.$TARGET.run.log
	java -Djava.io.tmpdir=$TMPDIR/$OUT.$TARGET -jar /$PICARD/picard.jar SortVcf \
		I=$OUTDIR/$OUT.$TARGET.filtered.snps.vcf.gz \
		I=$OUTDIR/$OUT.$TARGET.filtered.indels.vcf.gz \
		O=$OUTDIR/$OUT.$TARGET.filtered.sort.vcf.gz \
		TMP_DIR=$TMPDIR/$OUT.$TARGET
fi

# clean up
if [ -s $OUTDIR/$OUT.$TARGET.filtered.sort.vcf.gz.tbi ]; then
	echo -e "$(date)\nSorted VCFs found, cleaning up ...\n" &>> $LOGDIR/$OUT.$TARGET.run.log
	rm $OUTDIR/$OUT.$TARGET.filtered.snps.vcf.gz* &
	rm $OUTDIR/$OUT.$TARGET.filtered.indels.vcf.gz* &
	rm $OUTDIR/$OUT.$TARGET.indels.vcf.gz* &
	rm $OUTDIR/$OUT.$TARGET.snps.vcf.gz* &
	rm -r $TMPDIR/$OUT.$TARGET
fi

wait




