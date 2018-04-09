#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition Lewis  # for jobs < 2hrs try 'General'
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MPI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=20 # cores per task
#SBATCH --mem-per-cpu=8G  # memory per core (default is 1GB/core)
#SBATCH --time 2-00:00  # days-hours:minutes
#SBATCH --qos=normal
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
# IDln    lane no
# PL: instrument
# LB: library name
# SM: sample name
# PI: insert size
# R1: filename for read 1 (no .gz) 
# R2: filename for read 2 (no .gz)
# D1: dir for read 1 / dir name for uBAM
# D2: dir for read 2 / file name for uBAM
# REF: ref genome base
# REFDIR: dir of ref genome
# IDX: indexed ref base
# IDXDIR: dir of indexed ref genome
# FQDIR: dir for uncompressed fastq files
# MAPDIR: output dir
# LOGDIR: log file dir
# MEDIR: metrics file dir

# here we can print the RG info

#rsync -av --no-p --no-o --no-g -I --size-only lewis4-dtn.rnet.missouri.edu://storage/htc/lyonslab/washU_unprocessed/180119_84560625518588

# START=$(date +%g%m%d%H%M%S) sbatch --array=$(ls ./*sample_sheet.txt) map_libraries.sh

#-------------------------------------------------------------------------------

module load bwa/bwa-0.7.17
module load samtools/samtools-1.7
module load pigz/pigz-2.4
module load gatk/gatk-4.0.1.1

THREADS=20

SAMPLE_SHEET=$SLURM_ARRAY_TASK_ID

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
	read IDfl IDln PL LB SM PI R1 R2 D1 D2 REF REFDIR IDX IDXDIR FQDIR MAPDIR LOGDIR MEDIR <<< $(sed "${ROW}q;d" $SAMPLE_SHEET)

	# create a log dir so everything has somewhere to go
	if [ -d $LOGDIR/$START ]
        then	
                echo "$LOGDIR/$START found" &>> $LOGDIR/$START/$SM.run.log
        else
                mkdir -p $LOGDIR/$START
                touch $LOGDIR/$START/$SM.run.log
                echo "LOGDIR not found, making new one called $LOGDIR/$START" &>> $LOGDIR/$START/$SM.run.log
        fi


	# check if fastq files exist
	if [[ $D2 = *".bam" ]]; 
		then
			if [[ -e $D1/$D2 ]]; 
				then
					echo bam file $D2/$D1 found &>> $LOGDIR/$START/$SM.run.log
				else
					echo bam file $D2/$D1 not found, exiting &>> $LOGDIR/$START/$SM.run.log
					exit
				fi
		else
			if [[ -e $D1/$R1.gz && -e $D2/$R2.gz ]]; 
				then
					echo fastq files $D1/$R1.gz and $D2/$R2.gz found &>> $LOGDIR/$START/$SM.run.log
				else
					echo fastq files $D1/$R1.gz and $D2/$R2.gz found, exiting &>> $LOGDIR/$START/$SM.run.log
					exit
				fi

	# check if ref exists
	if [[ -e $REFDIR/$REF ]]; 
				then
					echo ref file $REFDIR/$REF found &>> $LOGDIR/$START/$SM.run.log
				else
					echo ref file $REFDIR/$REF not found, exiting &>> $LOGDIR/$START/$SM.run.log
					exit
				fi

	# check if index exists
	if [[ -e $IDXDIR/$IDX ]]; 
				then
					echo index file $IDXDIR/$IDX found &>> $LOGDIR/$START/$SM.run.log
				else
					echo index file $IDXDIR/$IDX not found, exiting &>> $LOGDIR/$START/$SM.run.log
					exit
				fi
	echo "\n\n"
done

echo "All index, ref and data files were found, starting RG alignment\n" &>> $LOGDIR/$START/$SM.run.log



# begin processing data
for ROW in $(seq 1 $(wc -l < $SAMPLE_SHEET))
do
	read IDfl IDln PL LB SM PI R1 R2 D1 D2 REF REFDIR IDX IDXDIR FQDIR MAPDIR LOGDIR MEDIR <<< $(sed "${ROW}q;d" $SAMPLE_SHEET)
	

    # state read group
	RG="@RG\tID:$IDfl:$IDln\tPL:$PL\tLB:$LB\tSM:$SM\tPI:$PI"
	echo -e "\n\n\n"At $(date) processing read group:"\n"$RG &>> $LOGDIR/$START/$SM.run.log


	if [ -d $FQDIR ]
		then
        	echo "FQDIR found" &>> $LOGDIR/$START/$SM.run.log
		else
			echo "FQDIR not found, making new one called $FQDIR" &>> $LOGDIR/$START/$SM.run.log
			mkdir $FQDIR	
		fi

	if [ -d $MAPDIR ]
        then
                echo "MAPDIR found" &>> $LOGDIR/$START/$SM.run.log
        else
                echo "MAPDIR not found, making new one called $MAPDIR" &>> $LOGDIR/$START/$SM.run.log
                mkdir $MAPDIR
        fi

	if [ -d $MEDIR ]
        then
                echo "MEDIR found" &>> $LOGDIR/$START/$SM.run.log
        else
                echo "MEDIR not found, making new one called $MEDIR" &>> $LOGDIR/$START/$SM.run.log
                mkdir $MEDIR
        fi

	
	# uncompress the files
	if [[ $D2 = *".bam" ]]; 
		then
			echo data is storred in unaligned bam format, converting to fastq &>> $LOGDIR/$START/$SM.run.log
			samtools fastq --threads $THREADS -1 $FQDIR/$R1 -2 $FQDIR/$R2 $D1/$D2 &>> $LOGDIR/$START/$SM.run.log
		else
			echo data is likely storred as compressed fastq, uncompressing &>> $LOGDIR/$START/$SM.run.log
			pigz -cd -p $THREADS $D1/$R1.gz > $FQDIR/$R1
			pigz -cd -p $THREADS $D2/$R2.gz > $FQDIR/$R2
		fi

	# check if file succesfully uncompressed
	if [[ -s $FQDIR/$R1 && -s $FQDIR/$R2 ]]; 
		then
			echo uncompressed read pair files found &>> $LOGDIR/$START/$SM.run.log
		else
			echo one or both uncompressed files are not found or are empty, exiting &>> $LOGDIR/$START/$SM.run.log
			exit
		fi


	echo begin mapping &>> $LOGDIR/$START/$SM.run.log
	# perform mapping
	(bwa mem -M -R $RG -t $THREADS $IDXDIR/$IDX $FQDIR/$R1 $FQDIR/$R2 | samtools view -Sb - > $MAPDIR/$SM.$LB.$IDfl.$IDln.bam) 2> $LOGDIR/$START/$SM.$LB.$IDfl.$IDln.aln.log

	#check for bam files and remove fastq files once reads are mapped
	if [ -s $MAPDIR/$SM.$LB.$IDfl.$IDln.bam ]
	then 
		echo bam file found, removing $FQDIR/{$R1,$R2} &>> $LOGDIR/$START/$SM.run.log
		rm -r $FQDIR/{$R1,$R2}
	else
		echo bam file not found or is empty, exiting &>> $LOGDIR/$START/$SM.run.log
		exit
	fi

	echo end mapping, begin sort &>> $LOGDIR/$START/$SM.run.log

	# sort and mark optical duplicates
	samtools sort --threads $THREADS -o $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam $MAPDIR/$SM.$LB.$IDfl.$IDln.bam &>> $LOGDIR/$START/$SM.run.log

	if [ -s $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam ]
	then
		echo sorted bam found, removing $MAPDIR/$SM.$LB.$IDfl.$IDln.bam &>> $LOGDIR/$START/$SM.run.log
		rm $MAPDIR/$SM.$LB.$IDfl.$IDln.bam
	else
		echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM.run.log
		exit
	fi

	echo end sort, begin duplicate marking &>> $LOGDIR/$START/$SM.run.log

	gatk MarkDuplicates -I=$MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam -O=$MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.markedDup.bam -M=$MEDIR/$SM.$LB.$IDfl.$IDln.metrics -OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 2> $LOGDIR/$START/$SM.$LB.$IDfl.$IDln.markDup.log
	# check for sorted bams and rm unsorted bams
	# using recomended distance for patterend flowcells, default is 100.

	if [ -s $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.markedDup.bam ]
        then
                echo duplicate marked bam found, removing $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam &>> $LOGDIR/$START/$SM.run.log
                rm $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam
        else
                echo duplicate marked bam not found or is empty, exiting &>> $LOGDIR/$START/$SM.run.log
                exit
        fi
    echo end duplicate marking &>> $LOGDIR/$START/$SM.run.log
done


echo -e "read groups have all been processed, begin individual sample processing\n" &>> $LOGDIR/$START/$SM.run.log
# here we merge files, this is done on sample ID, which is the first name ID

for SAMPLE in $(awk '{print $5}' $SAMPLE_SHEET | uniq)
do
	echo -e "\n" begin merge for sample $SAMPLE "\n" &>> $LOGDIR/$START/$SM.run.log
	samtools merge --threads $THREADS $MAPDIR/$SAMPLE.bam $MAPDIR/$SAMPLE.*.markedDup.bam &>> $LOGDIR/$START/$SM.run.log
	if [ -s $MAPDIR/$SAMPLE.bam ]
	then 
		echo merged bam found, removing non-merged bams &>> $LOGDIR/$START/$SM.run.log
		rm $MAPDIR/$SAMPLE.*.markedDup.bam
	else
		echo merged bam not found or is empty, exiting &>> $LOGDIR/$START/$SM.run.log
		exit
	fi
	

	echo -e "\n" end merge for sample $SAMPLE, begin sort for sample $SAMPLE "\n" &>> $LOGDIR/$START/$SM.run.log
	# sort and mark PCR duplicates
        samtools sort --threads $THREADS -o $MAPDIR/$SAMPLE.sorted.bam $MAPDIR/$SAMPLE.bam &>> $LOGDIR/$START/$SM.run.log

        if [ -s $MAPDIR/$SAMPLE.sorted.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SAMPLE.bam &>> $LOGDIR/$START/$SM.run.log
                rm $MAPDIR/$SAMPLE.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM.run.log
                exit
        fi
    
    echo -e "\n" end sort for sample $SAMPLE, begin duplicate marking for sample $SAMPLE "\n" &>> $LOGDIR/$START/$SM.run.log
        gatk MarkDuplicates -I=$MAPDIR/$SAMPLE.sorted.bam -O=$MAPDIR/$SAMPLE.sorted.markedDup.bam -M=$MEDIR/$SAMPLE.metrics 2> $LOGDIR/$START/$SAMPLE.markDup.log

        # check for sorted bams and rm unsorted bams
        if [ -s $MAPDIR/$SAMPLE.sorted.markedDup.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SAMPLE.sorted.bam &>> $LOGDIR/$START/$SM.run.log
                rm $MAPDIR/$SAMPLE.sorted.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/$SM.run.log
                exit
        fi
    echo -e "\n" end duplicate marking for sample $SAMPLE "\n" &>> $LOGDIR/$START/$SM.run.log
done

echo mapping completed for sample sheet: $SAMPLE_SHEET &>> $LOGDIR/$START/$SM.run.log
cat $SAMPLE_SHEET &>> $LOGDIR/$START/$SM.run.log
echo done
