#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p BioCompute
#SBATCH -J varcall_raw_bam
#SBATCH -o varcall_raw_bam.o%j
#SBATCH -e varcall_raw_bam.e%j
#SBATCH --mem 200G
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -t 2-00:00


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#----------------------------------------------------------

#----------------------------------------------------------
# CONFIG FOR JOB VAR
#----------------------------------------------------------
BAM_DIR=/home/buckleyrm/storage.lyonslab/big.test.preprocess/symlink.bam
GVCF_DIR=/home/buckleyrm/storage.lyonslab/big.test.preprocess/gvcf
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa
#----------------------------------------------------------

#invoke with 
## sbatch --array=0-$(( $(ls <bam dir> | grep bam$ | wc -l | cut -d " " -f 1 ) - 1)) call_var.sh

module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8
#module load picard-tools/picard-tools-2.1.1

# select cat to run on single node
LIST=($(ls $BAM_DIR | grep "bam$"))
CAT_BAM=${LIST[$SLURM_ARRAY_TASK_ID]}
#CAT_BAM=Fcat-22055-Chediak.sorted.markedDup.bam

echo ""
echo "Sample being analysed is $CAT_BAM"
echo ""

echo "user variables are bam_dir = $BAM_DIR, gvcf_dir $GVCF_DIR, ref = $REF"


java -Djava.io.tmpdir=$GVCF_DIR/tmp -XX:ParallelGCThreads=2 -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
-nct 26 \
-ERC GVCF \
-T HaplotypeCaller \
-R $REF \
-I $BAM_DIR/$CAT_BAM \
-o $GVCF_DIR/$CAT_BAM.g.vcf.gz

