#!/bin/bash
#SBATCH -p BioCompute,htc4,hpc5,Lewis
#SBATCH --account=general
#SBATCH -J Filter-exome
#SBATCH --mem 10G
#SBATCH -N1
#SBATCH -n2
#SBATCH -t 1-00:00
#SBATCH --output=Filter-exome.out

## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send

#------------------------------------------------------------------
#USER defined option
REF=/storage/htc/lyonslab/reference_files/fasta_genome/Felis_catus_9.0/Felis_catus_9.0.fa

GATK=/cluster/software/gatk/gatk-3.8/

INDIR=/storage/htc/warrenlab/users/alanarodney/exomes/gvcf_combine/out
IN=final_out
OUTDIR=/storage/htc/warrenlab/users/alanarodney/exomes/gvcf_combine/out/filtered_vcf
OUT=filtered_exome
TMP=/storage/htc/warrenlab/users/alanarodney/exomes/gvcf_combine/out/filtered_vcf/TMP

java -Djava.io.tmpdir=$TMP -Xmx9G -jar $GATK/GenomeAnalysisTK.jar \
                -T SelectVariants \
                -R $REF \
                -V $INDIR/$IN.vcf.gz \
                -selectType SNP \
                --log_to_file $TMP/$OUT.$TARGET.snp.log \
                -o $OUTDIR/$OUT.snp.vcf.gz

java -Djava.io.tmpdir=$TMP -Xmx9G -jar $GATK/GenomeAnalysisTK.jar \
                -T SelectVariants \
                -R $REF \
                -V $INDIR/$IN.vcf.gz \
                -selectType INDEL \
                --log_to_file $TMP/$OUT.$TARGET.indel.log \
                -o $OUTDIR/$OUT.indel.vcf.gz


java -Djava.io.tmpdir=$TMP -Xmx9G -jar $GATK/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $REF \
                -V $OUTDIR/$OUT.snp.vcf.gz \
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
                --log_to_file $TMP/$OUT.filtered.snp.log \
                -o $OUTDIR/$OUT.filtered.snp.vcf.gz

java -Djava.io.tmpdir=$TMP -Xmx9G -jar $GATK/GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $REF \
                -V $OUTDIR/$OUT.indel.vcf.gz \
                --filterExpression "QD < 2.0" \
                --filterName "QD" \
                --filterExpression "FS > 200.0" \
                --filterName "FS" \
                --filterExpression "SOR > 10.0" \
                --filterName "SOR" \
                --filterExpression "ReadPosRankSum < -20.0" \
                --filterName "ReadPosRankSum" \
                --log_to_file $TMP/$OUT.filtered.indel.log \
                -o $OUTDIR/$OUT.filtered.indel.vcf.gz

java -Djava.io.tmpdir=$TMP -Xmx9G -jar $PICARD/picard.jar MergeVcfs \
                I=$OUTDIR/$OUT.filtered.snp.vcf.gz \
                I=$OUTDIR/$OUT.filtered.indel.vcf.gz \
                O=$OUTDIR/$OUT.vcf.gz \
                TMP_DIR=$TMP

if [ -s $OUTDIR/$OUT.vcf.gz.tbi ]; then
        rm $OUTDIR/$OUT.filtered.snp.vcf.gz* &
        rm $OUTDIR/$OUT.filtered.indel.vcf.gz* &
        rm $OUTDIR/$OUT.indel.vcf.gz* &
        rm $OUTDIR/$OUT.snp.vcf.gz* &
fi


