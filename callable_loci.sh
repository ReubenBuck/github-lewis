#!/bin/bash
#----------------------------------------------------------
# CONFIG FOR SBATCH
#----------------------------------------------------------
#SBATCH -p BioCompute
#SBATCH --account=biocommunity
#SBATCH -J callable_loci
#SBATCH --output=Callable_loci-%A_%a.out
#SBATCH --mem 150G
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 2-00:00


## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#----------------------------------------------------------

#----------------------------------------------------------
# CONFIG FOR JOB VAR
#----------------------------------------------------------
BAM_DIR=/home/buckleyrm/storage.hpc.buckelyrm/printed_bams/Strict/recal/bams
OUT_DIR=/home/buckleyrm/storage.hpc.buckelyrm/callable
REF=/home/buckleyrm/storage.lyonslab/cat_ref/Felis_catus_9.0.fa
LOGDIR=/home/buckleyrm/storage.hpc.buckelyrm/callable/logs
TMPDIR=/storage/hpc/group/UMAG/WORKING/buckleyrm/TMP
#----------------------------------------------------------

#invoke with 
## sbatch --array=0-$(expr $(ls $BAM_DIR | grep ".bam$" | wc -l | cut -d " " -f 1) - 1) callable_loci.sh
module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8


# select cat to run on single node
LIST=($(ls -1 -X -r $BAM_DIR | grep "bam$"))
CAT_BAM=${LIST[$SLURM_ARRAY_TASK_ID]}

java -Djava.io.tmpdir=$TMPDIR/${CAT_BAM/.bam/} -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar \
        -T CallableLoci \
        -R $REF \
        -I $BAM_DIR/$CAT_BAM \
        -summary $OUT_DIR/${CAT_BAM/.bam/}.table.txt \
        -o $OUT_DIR/${CAT_BAM/.bam/}.callable_status.bed \
        --log_to_file $LOGDIR/${CAT_BAM/.bam/}.callable.log


