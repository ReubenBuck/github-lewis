#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p BioCompute,hpc5,hpc4
#SBATCH -J varcall_gatk_recal_bam
#SBATCH --output varcall_gatk_recal_bam.%A_%a.out
#SBATCH --mem 200G
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -t 7-00:00
#SBATCH --account=biocommunity
#SBATCH --qos=biolong


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#----------------------------------------------------------

#----------------------------------------------------------
# CONFIG FOR JOB VAR
#----------------------------------------------------------
BAM_DIR=/storage/htc/lyonslab/results/Felis_catus/bams
GVCF_DIR=/storage/htc/lyonslab/results/Felis_catus/gvcf
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa
LIST=/home/buckleyrm/storage.lyonslab/dom_cat_run/remaining_wilds/remaining_wilds.list
BAM_SUF="sort.markDup.realign.bam"
OUT_SUF=".g.vcf.gz"
#----------------------------------------------------------



#invoke with 
## sbatch --array=1-$(cat <cat.list> | wc -l | cut -d " " -f 1) call_var.sh

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8


sleep $((RANDOM % 10))


# select cat to run on single node
CAT_BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)

echo ""
echo "Sample being analysed is $CAT_BAM"
echo ""

echo "user variables are bam_dir = $BAM_DIR, gvcf_dir $GVCF_DIR, ref = $REF"

mkdir -p $GVCF_DIR/tmp/$CAT_BAM

java -Djava.io.tmpdir=$GVCF_DIR/tmp/$CAT_BAM -XX:ParallelGCThreads=2 -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nct 25 \
-ERC GVCF \
-T HaplotypeCaller \
-R $REF \
-I $BAM_DIR/$CAT_BAM/$CAT_BAM.$BAM_SUF \
-o $GVCF_DIR/$CAT_BAM.$OUT_SUF

