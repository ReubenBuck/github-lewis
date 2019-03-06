#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p BioCompute
#SBATCH --account=biocommunity
#SBATCH -J BQSR_Strict
#SBATCH --output=BQSR_Strict-%A_%a.out
#SBATCH --mem 150G
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 2-00:00


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#----------------------------------------------------------

#----------------------------------------------------------
# CONFIG FOR JOB VAR
#----------------------------------------------------------
BAM_DIR=/storage/htc/lyonslab/users/buckleyrm/BQSR/bam_sims
OUT_DIR=/storage/htc/lyonslab/users/buckleyrm/BQSR/tables/Strict
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa
DATABASE=/home/buckleyrm/storage.lyonslab/v9_vcf/Strict/select_sites_Strict.vcf.gz
BAMOUTDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/printed_bams/Strict
LOGDIR=/storage/htc/lyonslab/users/buckleyrm/BQSR/logs/Strict
TMPDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/TMP
#----------------------------------------------------------

#invoke with 
## sbatch --array=0-$(expr $(ls $BAM_DIR | grep ".bam$" | wc -l | cut -d " " -f 1) - 1) BQSR.sh
module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
module load samtools/samtools-1.9 
module load r/r-3.4.2


# select cat to run on single node
LIST=($(ls -1 -X -r $BAM_DIR | grep "bam$"))
CAT_BAM=${LIST[$SLURM_ARRAY_TASK_ID]}


echo ""
echo "Sample being recalibrated is $CAT_BAM"
echo ""

echo "user variables are bam_dir = $BAM_DIR, out_dir $OUT_DIR, ref = $REF, known sites = $DATABASE"

sleep $((RANDOM % 10))

if [ -s $BAMOUTDIR/$CAT_BAM.bai ]; then
	echo header has been updated and indexed
else
	echo updating header and indexing file
	samtools view -H $BAM_DIR/$CAT_BAM | sed -e "s/\tSM:${CAT_BAM/.bam/}/\tSM:${CAT_BAM/.bam/}\tPL:ILLUMINA/g" | samtools reheader -P - $BAM_DIR/$CAT_BAM > $BAMOUTDIR/$CAT_BAM
	samtools index -@ 20 $BAMOUTDIR/$CAT_BAM
fi

sleep $((RANDOM % 10))

if [ ! -d $TMPDIR/${CAT_BAM/.bam/} ]; then
	mkdir $TMPDIR/${CAT_BAM/.bam/}
fi

# first pass

if [ -s $OUT_DIR/$CAT_BAM.recal_data.table ]; then
	echo found first pass recalibration table
else
	echo creating first pass recalibration table
	java -Djava.io.tmpdir=$TMPDIR/${CAT_BAM/.bam/} -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
	-nct 20 \
	-T BaseRecalibrator \
	-R $REF \
	-I $BAMOUTDIR/$CAT_BAM \
	-knownSites $DATABASE \
	-o $OUT_DIR/$CAT_BAM.recal_data.table \
	--log_to_file $LOGDIR/$CAT_BAM.first.log
fi

# second pass
if [ -s $OUT_DIR/$CAT_BAM.post_recal_data.table ]; then
        echo found second pass recalibration table
else
        echo creating second pass recalibration table
	java -Djava.io.tmpdir=$TMPDIR/${CAT_BAM/.bam/} -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
	-nct 20 \
	-T BaseRecalibrator \
	-R $REF \
	-I $BAMOUTDIR/$CAT_BAM \
	-knownSites $DATABASE \
	-BQSR $OUT_DIR/$CAT_BAM.recal_data.table \
	-o $OUT_DIR/$CAT_BAM.post_recal_data.table \
	--log_to_file $LOGDIR/$CAT_BAM.second.log
fi

# Generate before and afer plots

if [ -s $BAMOUTDIR/${CAT_BAM/bam/recal.bam} ]; then
	echo found recalibrated bam
else
	echo printing reads for recalibrated bam
	java -Djava.io.tmpdir=$TMPDIR/${CAT_BAM/.bam/} -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
	-nct 20 \
	-T PrintReads \
	-R $REF \
	-I  $BAMOUTDIR/$CAT_BAM \
	-BQSR $OUT_DIR/$CAT_BAM.recal_data.table \
	-o $BAMOUTDIR/${CAT_BAM/bam/recal.bam} \
	--log_to_file $LOGDIR/$CAT_BAM.print_reads.log
fi


