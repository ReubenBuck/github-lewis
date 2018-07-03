#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p General
#SBATCH -J BQSR
#SBATCH -o BQSR.o%j
#SBATCH -e BQSR.e%j
#SBATCH --mem 2G
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 0-00:10


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#----------------------------------------------------------

#----------------------------------------------------------
# CONFIG FOR JOB VAR
#----------------------------------------------------------
BAM_DIR=/home/buckleyrm/storage.lyonslab/users/yuyos/FSEC/bams
OUT_DIR=/home/buckleyrm/storage.lyonslab/users/yuyos/FSEC/recal
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa
DATABASE=/home/buckleyrm/storage.lyonslab/v9_vcf/Unrelated_cats.Genotypefiltered.recode.vcf.gz
#----------------------------------------------------------

#invoke with 
## sbatch --array=1-$(ls $BAM_DIR | grep ".bam$" | wc -l | cut -d " " -f 1) BQSR.sh

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
#module load picard-tools/picard-tools-2.1.1

# select cat to run on single node
LIST=($(ls -1 -X -r $BAM_DIR | grep "bam$"))
CAT_BAM=${LIST[$SLURM_ARRAY_TASK_ID]}


echo ""
echo "Sample being recalibrated is $CAT_BAM"
echo ""

echo "user variables are bam_dir = $BAM_DIR, out_dir $OUT_DIR, ref = $REF, known sites = $DATABASE"

# first pass
java -Djava.io.tmpdir=/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nct 20 \
-T BaseRecalibrator \
-R $REF \
-I $BAM_DIR/$CAT_BAM \
-knownSites $DATABASE \
-o $OUT_DIR/$CAT_BAM.recal_data.table

# second pass
java -Djava.io.tmpdir=/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nct 20 \
-T BaseRecalibrator \
-R $REF \
-I $BAM_DIR/$CAT_BAM \
-knownSites $DATABASE \
-BQSR $OUT_DIR/$CAT_BAM.recal_data.table \
-o $OUT_DIR/$CAT_BAM.post_recal_data.table

# Generate before and afer plots
java -Djava.io.tmpdir=/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar
-nct 20 \
-T AnalyzeCovariates \
-R $REF \
-before $OUT_DIR/$CAT_BAM.recal_data.table \
-after $OUT_DIR/$CAT_BAM.post_recal_data.table \
-plots $OUT_DIR/$CAT_BAM.recal_plot.pdf

java -Djava.io.tmpdir=/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar
-nct 20 \
-T PrintReafd \
-R $REF \
-I  $BAM_DIR/$CAT_BAM \
-BQSR $CAT_BAM.recal_data.table \
-o $BAM_DIR/${CAT_BAM/bam/recal.bam}





