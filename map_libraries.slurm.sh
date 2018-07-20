#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition BioCompute  # for jobs < 2hrs try 'General'
#SBATCH -N1
#SBATCH -n52 # cores 
#SBATCH --mem 450G  # memory 
#SBATCH -t 2-00:00  # days-hours:minutes
#SBATCH --account=biocommunity  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=map_lib
#SBATCH --output=map_results-%A_%a.out  # %j is the unique jobID
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

# START=$(date +%g-%m-%d_%H.%M.%S) LIST=<samp_sheet_paths> sbatch --array=1-$(wc -l <samp_sheet_paths> | cut -d " " -f 1) map_libraries.slurm.sh

# IMPORTANT: REF directory requries a dir called target_loci. Inside this dir there needs to be files that contain intervals that cover the entire genome.
# The way the intervals are split across the files in "target_loci" will dictate how the data is split up for scatter gather approach in indel realignment.

#-------------------------------------------------------------------------------

module load bwa/bwa-0.7.17
module load samtools/samtools-1.7
module load pigz/pigz-2.4
module load fastqc
module load picard-tools/picard-tools-2.1.1
module load java/openjdk/java-1.8.0-openjdk
module load gatk/gatk-3.8

THREADS=52


echo $SLURM_ARRAY_TASK_ID

SAMPLE_SHEET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $LIST)

echo $SAMPLE_SHEET


if [ -z ${START+x} ];
	then 
		echo "START time is empty, check cmd line, exiting"
		exit
	else 
		echo "START time is set to $START, process is begining"
	fi


# before we start cheking everything we set up a wait time
sleep $((RANDOM % 10))
 

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


	if [ -d $FQDIR/$SM ]
		then
        	echo "FQDIR found" &>> $LOGDIR/$START/$SM/$SM.run.log
		else
			echo "FQDIR not found, making new one called $FQDIR" &>> $LOGDIR/$START/$SM/$SM.run.log
			mkdir -p $FQDIR/$SM/tmp
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
                echo "MEDIR not found, making new one called $MEDIR/$SM" with tmp dir &>> $LOGDIR/$START/$SM/$SM.run.log
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
			samtools fastq --threads $THREADS -1 $FQDIR/$SM/$R1 -2 $FQDIR/$SM/$R2 $D1/$D2 &>> $LOGDIR/$START/$SM/$SM.run.log

			#test params
			#samtools view -h $D1/$D2 | head -n 15000 | samtools view -Sb - > $FQDIR/$SM/tmp.$D2
			#samtools fastq --threads $THREADS -1 $FQDIR/$SM/$R1 -2 $FQDIR/$SM/$R2 $FQDIR/$SM/tmp.$D2 &>> $LOGDIR/$START/$SM/$SM.run.log
			
		else
			echo data is likely storred as compressed fastq, uncompressing &>> $LOGDIR/$START/$SM/$SM.run.log
			pigz -cd -p $THREADS $D1/$R1.gz > $FQDIR/$SM/$R1
			pigz -cd -p $THREADS $D2/$R2.gz > $FQDIR/$SM/$R2
			
			#test params
			#pigz -cd -p $THREADS $D1/$R1.gz | head -n 40000 > $FQDIR/$SM/$R1
                        #pigz -cd -p $THREADS $D2/$R2.gz | head -n 40000 > $FQDIR/$SM/$R2
		fi

	# check if file succesfully uncompressed
	if [[ -s $FQDIR/$SM/$R1 && -s $FQDIR/$SM/$R2 ]]; 
		then
			echo uncompressed read pair files found &>> $LOGDIR/$START/$SM/$SM.run.log
		else
			echo one or both uncompressed files are not found or are empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
			exit
		fi

	# check quality with fastqc
	echo running fastqc quality checks ... &>> $LOGDIR/$START/$SM/$SM.run.log
	fastqc -o $QUALDIR/$SM $FQDIR/$SM/$R1 $FQDIR/$SM/$R2

	echo begin mapping &>> $LOGDIR/$START/$SM/$SM.run.log
	# perform mapping
	(bwa mem -M -R $RG -t $THREADS $IDX $FQDIR/$SM/$R1 $FQDIR/$SM/$R2 | samtools view -Sb - > $MAPDIR/$SM/$SM.$ROW.bam) 2> $LOGDIR/$START/$SM/$SM.$ROW.aln.log

	#check for bam files and remove fastq files once reads are mapped
	if [ -s $MAPDIR/$SM/$SM.$ROW.bam ]
	then 
		echo bam file found, removing $FQDIR/$SM/{$R1,$R2} &>> $LOGDIR/$START/$SM/$SM.run.log
		rm -r $FQDIR/$SM/{$R1,$R2}
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
		rm $(eval echo $MAPDIR/$SM/$SM.{1..$(wc -l $SAMPLE_SHEET | cut -d " " -f 1)}.bam)
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
                rm $MAPDIR/$SM/$SM.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                exit
        fi
    
    echo -e "\n" end sort for sample $SM, begin duplicate marking for sample $SM "\n" &>> $LOGDIR/$START/$SM/$SM.run.log

java -Djava.io.tmpdir=$FQDIR/$SM/tmp -jar /cluster/software/picard-tools/picard-tools-2.1.1/picard.jar MarkDuplicates TMP_DIR=$FQDIR/$SM/tmp INPUT=$MAPDIR/$SM/$SM.sorted.bam OUTPUT=$MAPDIR/$SM/$SM.sorted.markedDup.bam METRICS_FILE=$MEDIR/$SM/$SM.metrics OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2> $LOGDIR/$START/$SM/$SM.markDup.log

        # check for sorted bams and rm unsorted bams
        if [ -s $MAPDIR/$SM/$SM.sorted.markedDup.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SM/$SM.sorted.bam &>> $LOGDIR/$START/$SM/$SM.run.log
                rm $MAPDIR/$SM/$SM.sorted.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                exit
        fi
    echo -e "\n" end duplicate marking for sample $SM "\n" &>> $LOGDIR/$START/$SM/$SM.run.log

samtools index -@ $THREADS $MAPDIR/$SM/$SM.sorted.markedDup.bam

samtools stats --threads $THREADS $MAPDIR/$SM/$SM.sorted.markedDup.bam > $QUALDIR/$SM/$SM.sorted.markedDup.bam.bc
mkdir $QUALDIR/$SM/markDupQualPlot
plot-bamstats -p $QUALDIR/$SM/markDupQualPlot/$SM.markDup $QUALDIR/$SM/$SM.sorted.markedDup.bam.bc

echo mapping completed for sample sheet: $SAMPLE_SHEET &>> $LOGDIR/$START/$SM/$SM.run.log

echo begin indel realignment, generate targets &>> $LOGDIR/$START/$SM/$SM.run.log

java -Djava.io.tmpdir=$FQDIR/$SM/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar -nt $THREADS -T RealignerTargetCreator -R $REF -I $MAPDIR/$SM/$SM.sorted.markedDup.bam -o $MEDIR/$SM/$SM.indelTarget.intervals 2> $LOGDIR/$START/$SM/$SM.indelTargetCreator.log

echo targets generated, begin realignment &>> $LOGDIR/$START/$SM/$SM.run.log

# run indel realignment for avialbale threads only

N=$(expr $THREADS - 1)
for TARGET in $(ls ${REF%/*}/target_loci/); do
	(
	sleep $((RANDOM % 20))
	echo begin realignment for ${TARGET%\.intervals} at $(date) &>> $LOGDIR/$START/$SM/$SM.run.log
	java -Djava.io.tmpdir=$FQDIR/$SM/tmp -jar /cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $MAPDIR/$SM/$SM.sorted.markedDup.bam -targetIntervals $MEDIR/$SM/$SM.indelTarget.intervals -L ${REF%/*}/target_loci/$TARGET -o $MAPDIR/$SM/$SM.${TARGET%\.intervals}.sorted.markedDup.realigned.bam 2> $LOGDIR/$START/$SM/$SM.${TARGET%\.intervals}.realign.log 
	) &

	# will make sure that we do not open up more threads then there are available
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
		wait
	fi

done

# wait until all indels have been realigned
wait

# check if all target files got realigned

for TARGET in $(ls ${REF%/*}/target_loci/); do
        if [ -s $MAPDIR/$SM/$SM.${TARGET%\.intervals}.sorted.markedDup.realigned.bam ]
                then
                        echo realigned bam for ${TARGET%\.intervals} found, continuing &>> $LOGDIR/$START/$SM/$SM.run.log
                else
                        echo realigned bam for ${TARGET%\.intervals} not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                        exit
                fi
done

echo '\nall targets were realigned for $SM at $(date), removing non relaigned bams and concatenating\n' &>> $LOGDIR/$START/$SM/$SM.run.log

rm $MAPDIR/$SM/$SM.sorted.markedDup.bam $MAPDIR/$SM/$SM.sorted.markedDup.bam.bai

# is it a good idea to merge at this stage?, yes since first pass of score reclibration requires whole genome
samtools cat -o $MAPDIR/$SM/$SM.realign.bam $(ls ${REF%/*}/target_loci/ | sed "s|^|$MAPDIR/$SM/$SM.|g" | sed "s/.intervals/.sorted.markedDup.realigned.bam/g")

# remove all of the old files
if [ -s $MAPDIR/$SM/$SM.realign.bam ]
	then
        	echo concatenated realigned bam for $SM found, remove target files &>> $LOGDIR/$START/$SM/$SM.run.log
		rm $(ls ${REF%/*}/target_loci/ | sed "s|^|$MAPDIR/$SM/$SM.|g" | sed "s/.intervals/.sorted.markedDup.realigned.bam/g")
		rm $(ls ${REF%/*}/target_loci/ | sed "s|^|$MAPDIR/$SM/$SM.|g" | sed "s/.intervals/.sorted.markedDup.realigned.bai/g")
        else
        	echo concatenated realigned bam for $SM not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                exit
        fi



# then probably sort and index the single file


        echo -e '\nend concatenate for sample $SM, begin sort for sample $SM \n' &>> $LOGDIR/$START/$SM/$SM.run.log
        # sort and mark PCR duplicates
        samtools sort --threads $THREADS -o $MAPDIR/$SM/$SM.sort.markDup.realign.bam $MAPDIR/$SM/$SM.realign.bam &>> $LOGDIR/$START/$SM/$SM.run.log

        if [ -s $MAPDIR/$SM/$SM.sort.markDup.realign.bam ]
        then
                echo sorted realigned bam found, removing $MAPDIR/$SM/$SM.realign.bam &>> $LOGDIR/$START/$SM/$SM.run.log
                rm $MAPDIR/$SM/$SM.realign.bam
        else
                echo sorted realigned bam not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
                exit
        fi

samtools index -@ $THREADS $MAPDIR/$SM/$SM.sort.markDup.realign.bam

samtools stats --threads $THREADS $MAPDIR/$SM/$SM.sort.markDup.realign.bam > $QUALDIR/$SM/$SM.sort.markDup.realign.bam.bc
mkdir $QUALDIR/$SM/realignQualPlot
plot-bamstats -p $QUALDIR/$SM/realignQualPlot/$SM.realign $QUALDIR/$SM/$SM.sort.markDup.realign.bam.bc

cat $SAMPLE_SHEET &>> $LOGDIR/$START/$SM/$SM.run.log


echo done
