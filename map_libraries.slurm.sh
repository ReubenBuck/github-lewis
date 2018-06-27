#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition Interactive  # for jobs < 2hrs try 'General'
#SBATCH -N1
#SBATCH -n3 # cores 
#SBATCH --mem 10G  # memory 
#SBATCH -t 0-00:20  # days-hours:minutes
#SBATCH --account=general  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=map_lib
#SBATCH --output=results-%j.out  # %j is the unique jobID
#
## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#
## array options
#
#-------------------------------------------------------------------------------

### NOTES


# a smaple sheet needs to be fitted to an individual sample, this way we can feed the sheets into the pipeline and set up an array
# this script takes care of the first part of preprocessing.
# the next few stages can be a little more complex
# ie base recalibration

# this will be useful for getting soem of the required samples into bam files. 

# input is a tab delimited file with the following fields:
# IDfl: flowcell ID
# IDln: lane no
# LB: library name
# SM: sample name
# R1: filename for read 1 (no .gz) 
# R2: filename for read 2 (no .gz)
# D1: dir for read 1 / dir name for uBAM
# D2: dir for read 2 / file name for uBAM
# REF: ref genome base
# IDX: indexed ref base
# FQDIR: dir for uncompressed fastq files
# MAPDIR: output dir
# LOGDIR: log file dir
# MEDIR: metrics file dir
# QUALDIR: directory for quality metrics

# START=$(date +%g%m%d%H%M%S) LIST=<samp_sheet_paths> sbatch --array=1-$(wc -l <samp_sheet_paths> | cut -d " " -f 1) map_libraries.slurm.sh
#-------------------------------------------------------------------------------

module load bwa/bwa-0.7.17
module load samtools/samtools-1.7
module load pigz/pigz-2.4
module load fastqc
module load picard-tools/picard-tools-2.1.1
module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

THREADS=3


echo $SLURM_ARRAY_TASK_ID

SAMPLE_SHEET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)

echo $SAMPLE_SHEET

#START=$(date +%g%m%d%H%M%S)

if [ -z ${START+x} ];
	then 
		echo "START time is empty, check cmd line, exiting"
		exit
	else 
		echo "START time is set to $START, process is begining"
	fi


# before we start cheking everything we set up a wait time
sleep $((RANDOM % 10))


# need to check if start was recorded 

# check if all files exist first. 
# ref files, idx files, fastq/bam files
for ROW in $(seq 1 $(wc -l < $SAMPLE_SHEET))
do
	read IDfl IDln LB SM R1 R2 D1 D2 REF IDX FQDIR MAPDIR LOGDIR MEDIR QUALDIR <<< $(sed "${ROW}q;d" $SAMPLE_SHEET)

	# create a log dir so everything has somewhere to go
	if [ -d $LOGDIR/$START/$SM ]
        then	
                echo "$LOGDIR/$START/$SM found" &>> $LOGDIR/$START/$SM/$SM.run.log
        else
                mkdir -p $LOGDIR/$START/$SM
                touch $LOGDIR/$START/$SM/$SM.run.log
                echo "LOGDIR not found, making new one called $LOGDIR/$START/$SM" &>> $LOGDIR/$START/$SM/$SM.run.log
        fi


	# check if fastq files exist
	if [[ $D2 = *".bam" ]]; 
		then
			if [[ -e $D1/$D2 ]]; 
				then
					echo bam file $D2/$D1 found &>> $LOGDIR/$START/$SM/$SM.run.log
				else
					echo bam file $D2/$D1 not found, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
					exit
				fi
		else
			if [[ -e $D1/$R1.gz && -e $D2/$R2.gz ]]; 
				then
					echo fastq files $D1/$R1.gz and $D2/$R2.gz found &>> $LOGDIR/$START/$SM/$SM.run.log
				else
					echo fastq files $D1/$R1.gz and $D2/$R2.gz found, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
					exit
				fi
	fi

	# check if ref exists
	if [[ -e $REF ]]; 
				then
					echo ref file $REF found &>> $LOGDIR/$START/$SM/$SM.run.log
				else
					echo ref file $REF not found, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
					exit
				fi

	# check if index exists
	if [[ -e $IDX ]]; 
				then
					echo index file $IDX found &>> $LOGDIR/$START/$SM/$SM.run.log
				else
					echo index file $IDX not found, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
					exit
				fi

	echo '\n\n' &>> $LOGDIR/$START/$SM/$SM.run.log
done

echo "All index, ref and data files were found, starting RG alignment\n" &>> $LOGDIR/$START/$SM/$SM.run.log



# begin processing data
for ROW in $(seq 1 $(wc -l < $SAMPLE_SHEET))
do
	read IDfl IDln LB SM R1 R2 D1 D2 REF IDX FQDIR MAPDIR LOGDIR MEDIR QUALDIR <<< $(sed "${ROW}q;d" $SAMPLE_SHEET)
	

    # state read group
	RG="@RG\tID:$LB:$IDfl:$IDln\tLB:$LB\tSM:$SM"
	echo -e "\n\n\n"At $(date) processing read group:"\n"$RG &>> $LOGDIR/$START/$SM/$SM.run.log


	if [ -d $FQDIR ]
		then
        	echo "FQDIR found" &>> $LOGDIR/$START/$SM/$SM.run.log
		else
			echo "FQDIR not found, making new one called $FQDIR" &>> $LOGDIR/$START/$SM/$SM.run.log
			mkdir $FQDIR	
		fi

	if [ -d $MAPDIR/$SM ]
        then
                echo "MAPDIR found" &>> $LOGDIR/$START/$SM/$SM.run.log
        else
                echo "MAPDIR not found, making new one called $MAPDIR/$SM" &>> $LOGDIR/$START/$SM/$SM.run.log
                mkdir -p $MAPDIR/$SM
        fi

	if [ -d $MEDIR/$SM ]
        then
                echo "MEDIR found" &>> $LOGDIR/$START/$SM/$SM.run.log
        else
                echo "MEDIR not found, making new one called $MEDIR/$SM" &>> $LOGDIR/$START/$SM/$SM.run.log
                mkdir -p $MEDIR/$SM
        fi

    if [ -d $QUALDIR/$SM ]
        then
                echo "QUALDIR found" &>> $LOGDIR/$START/$SM/$SM.run.log
        else
                echo "QUALDIR not found, making new one called $QUALDIR/$SM" &>> $LOGDIR/$START/$SM/$SM.run.log
                mkdir -p $QUALDIR/$SM
        fi

	
	# uncompress the files
	if [[ $D2 = *".bam" ]]; 
		then
			echo data is storred in unaligned bam format, converting to fastq &>> $LOGDIR/$START/$SM/$SM.run.log
			samtools view -h $D1/$D2 | head -n 105000 | samtools view -Sb - > $FQDIR/tmp.$D2
			samtools fastq --threads $THREADS -1 $FQDIR/$R1 -2 $FQDIR/$R2 $FQDIR/tmp.$D2 &>> $LOGDIR/$START/$SM/$SM.run.log
			rm $FQDIR/tmp.$D2
			#samtools fastq --threads $THREADS -1 $FQDIR/$R1 -2 $FQDIR/$R2 $D1/$D2 &>> $LOGDIR/$START/$SM/$SM.run.log
		else
			echo data is likely storred as compressed fastq, uncompressing &>> $LOGDIR/$START/$SM/$SM.run.log
			pigz -cd -p $THREADS $D1/$R1.gz | head -n 400000 > $FQDIR/$R1
			pigz -cd -p $THREADS $D2/$R2.gz | head -n 400000 > $FQDIR/$R2
		fi

	# check if file succesfully uncompressed
	if [[ -s $FQDIR/$R1 && -s $FQDIR/$R2 ]]; 
		then
			echo uncompressed read pair files found &>> $LOGDIR/$START/$SM/$SM.run.log
		else
			echo one or both uncompressed files are not found or are empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
			exit
		fi

	# check quality with fastqc
	echo running fastqc quality checks ... &>> $LOGDIR/$START/$SM/$SM.run.log
	#fastqc -o $QUALDIR/$SM $FQDIR/$R1 $FQDIR/$R2

	echo begin mapping &>> $LOGDIR/$START/$SM/$SM.run.log
	# perform mapping
	(bwa mem -M -R $RG -t $THREADS $IDX $FQDIR/$R1 $FQDIR/$R2 | samtools view -Sb - > $MAPDIR/$SM/$SM.$ROW.bam) 2> $LOGDIR/$START/$SM/$SM.$ROW.aln.log

	#check for bam files and remove fastq files once reads are mapped
	if [ -s $MAPDIR/$SM/$SM.$ROW.bam ]
	then 
		echo bam file found, removing $FQDIR/{$R1,$R2} &>> $LOGDIR/$START/$SM/$SM.run.log
		#rm -r $FQDIR/{$R1,$R2}
	else
		echo bam file not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
		exit
	fi

	echo end mapping &>> $LOGDIR/$START/$SM/$SM.run.log

done


echo -e "read groups have all been processed, begin individual sample processing\n" &>> $LOGDIR/$START/$SM/$SM.run.log
# here we merge files, this is done on sample ID, which is the first name ID


	echo -e "\n" begin merge for sample $SM "\n" &>> $LOGDIR/$START/$SM/$SM.run.log

	samtools merge -c -f --threads $THREADS $MAPDIR/$SM/$SM.bam $(eval echo $MAPDIR/$SM/$SM.{1..$(wc -l $SAMPLE_SHEET | cut -d " " -f 1)}.bam) &>> $LOGDIR/$START/$SM/$SM.run.log
	if [ -s $MAPDIR/$SM/$SM.bam ]
	then 
		echo merged bam found, removing non-merged bams &>> $LOGDIR/$START/$SM/$SM.run.log
	#	rm $(eval echo $MAPDIR/$SM/$SM.{1..$(wc -l $SAMPLE_SHEET | cut -d " " -f 1)}.bam)
	else
		echo merged bam not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
		exit
	fi
	

	echo -e "\n" end merge for sample $SM, begin sort for sample $SM "\n" &>> $LOGDIR/$START/$SM/$SM.run.log
	# sort and mark PCR duplicates
        samtools sort --threads $THREADS -o $MAPDIR/$SM/$SM.sorted.bam $MAPDIR/$SM/$SM.bam &>> $LOGDIR/$START/$SM/$SM.run.log

        if [ -s $MAPDIR/$SM/$SM.sorted.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SM/$SM.bam &>> $LOGDIR/$START/$SM/$SM.run.log
         #       rm $MAPDIR/$SM/$SM.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                exit
        fi
    
    echo -e "\n" end sort for sample $SM, begin duplicate marking for sample $SM "\n" &>> $LOGDIR/$START/$SM/$SM.run.log

java -jar /cluster/software/picard-tools/picard-tools-2.1.1/picard.jar MarkDuplicates INPUT=$MAPDIR/$SM/$SM.sorted.bam OUTPUT=$MAPDIR/$SM/$SM.sorted.markedDup.bam METRICS_FILE=$MEDIR/$SM/$SM.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2> $LOGDIR/$START/$SM/$SM.markDup.log

        # check for sorted bams and rm unsorted bams
        if [ -s $MAPDIR/$SM/$SM.sorted.markedDup.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SM/$SM.sorted.bam &>> $LOGDIR/$START/$SM/$SM.run.log
         #       rm $MAPDIR/$SM/$SM.sorted.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                exit
        fi
    echo -e "\n" end duplicate marking for sample $SM "\n" &>> $LOGDIR/$START/$SM/$SM.run.log

samtools index -@ $THREADS $MAPDIR/$SM/$SM.sorted.markedDup.bam

samtools stats --threads $THREADS $MAPDIR/$SM/$SM.sorted.markedDup.bam > $QUALDIR/$SM/$SM.sorted.markedDup.bam.bc
plot-bamstats -p $QUALDIR/$SM/ $QUALDIR/$SM/$SM.sorted.markedDup.bam.bc

echo mapping completed for sample sheet: $SAMPLE_SHEET &>> $LOGDIR/$START/$SM/$SM.run.log

echo begin indel realignment, generate targets &>> $LOGDIR/$START/$SM/$SM.run.log

java -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar -nt $THREADS -T RealignerTargetCreator -R $REF -I $MAPDIR/$SM/$SM.sorted.markedDup.bam -o $MEDIR/$SM/$SM.indelTarget.intervals 2> $LOGDIR/$START/$SM/$SM.indelTargetCreator.log

echo targets generated, begin realignment &>> $LOGDIR/$START/$SM/$SM.run.log


# run indel realignment for avialbale threads only

N=$(expr $THREADS - 1)

for TARGET in $(ls ${REF%/*}/target_loci/); do
	
	(
	echo relignment of ${TARGET%\.intervals} begining at $(date) &>> $LOGDIR/$START/$SM/$SM.run.log
	java -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $MAPDIR/$SM/$SM.sorted.markedDup.bam -targetIntervals $MEDIR/$SM/$SM.indelTarget.intervals -L ${REF%/*}/target_loci/$TARGET -o $MAPDIR/$SM/$SM.${TARGET%\.intervals}.sorted.markedDup.realigned.bam 2> $LOGDIR/$START/$SM/$SM.${TARGET%\.intervals}.realign.log 
	) &

	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
		wait
	fi

done

wait




cat $SAMPLE_SHEET &>> $LOGDIR/$START/$SM/$SM.run.log

echo done
