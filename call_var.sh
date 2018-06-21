#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p Lewis
#SBATCH -J varcall_raw_bam
#SBATCH -o varcall_raw_bam.o%j
#SBATCH -e varcall_raw_bam.e%j
#SBATCH --mem 200G
#SBATCH -N 1
#SBATCH -n 22
#SBATCH -t 2-00:00


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#----------------------------------------------------------

#----------------------------------------------------------
# CONFIG FOR JOB VAR
#----------------------------------------------------------
BAM_DIR=/home/buckleyrm/storage.lyonslab/users/yuyos/FSEC/bams
GVCF_DIR=/home/buckleyrm/storage.lyonslab/users/yuyos/FSEC/calls
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa
#----------------------------------------------------------

#invoke with 
## BAM_DIR=/home/buckleyrm/storage.lyonslab/users/yuyos/FSEC/bams sbatch --array=$(ls $BAM_DIR) call_var.sh

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
#module load picard-tools/picard-tools-2.1.1

# select cat to run on single node
LIST=($(ls -1 -X -r $BAM_DIR))
CAT_BAM=${LIST[$SLURM_ARRAY_TASK_ID]}
#CAT_BAM=Fcat-22055-Chediak.sorted.markedDup.bam

echo ""
echo "Sample being analysed is $CAT_BAM"
echo ""

echo "user variables are bam_dir = $BAM_DIR, gvcf_dir $GVCF_DIR, ref = $REF"


java -Djava.io.tmpdir=/tmp -XX:ParallelGCThreads=2 -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nct 20 \
-ERC GVCF \
-T HaplotypeCaller \
-R $REF \
-I $BAM_DIR/$CAT_BAM \
--pcr_indel_model NONE \
-o $GVCF_DIR/$CAT_BAM.g.vcf.gz

