#!/bin/bash/

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

# SAMPLE_SHEET=<sample sheet> THREADS=<number of threads for mapping> bash map_libraries.sh


START=$(date +%g%m%d%H%M%S)

for ROW in $(seq 1 $(wc -l < $SAMPLE_SHEET))
#for ROW in $(seq 1 2)
do
	read IDfl IDln PL LB SM PI R1 R2 D1 D2 REF REFDIR IDX IDXDIR FQDIR MAPDIR LOGDIR MEDIR<<< $(sed "${ROW}q;d" $SAMPLE_SHEET)
	
	if [ -d $LOGDIR/$START ]
        then	
        		echo "begin next read group\n" &>> $LOGDIR/$START/run.log
                echo "$LOGDIR/$START found" &>> $LOGDIR/$START/run.log
        else
                mkdir -p $LOGDIR/$START
                touch $LOGDIR/$START/run.log
                echo "LOGDIR not found, making new one called $LOGDIR/$START" &>> $LOGDIR/$START/run.log
        fi

    # state read group
	RG="@RG\tID:$IDfl:$IDln\tPL:$PL\tLB:$LB\tSM:$SM\tPI:$PI"
	echo -e At $(date) processing read group:"\n"$RG &>> $LOGDIR/$START/run.log


	if [ -d $FQDIR ]
		then
        	echo "FQDIR found" &>> $LOGDIR/$START/run.log
		else
			echo "FQDIR not found, making new one called $FQDIR" &>> $LOGDIR/$START/run.log
			mkdir $FQDIR	
		fi

	if [ -d $MAPDIR ]
        then
                echo "MAPDIR found" &>> $LOGDIR/$START/run.log
        else
                echo "MAPDIR not found, making new one called $MAPDIR" &>> $LOGDIR/$START/run.log
                mkdir $MAPDIR
        fi

	if [ -d $MEDIR ]
        then
                echo "MEDIR found" &>> $LOGDIR/$START/run.log
        else
                echo "MEDIR not found, making new one called $MEDIR" &>> $LOGDIR/$START/run.log
                mkdir $MEDIR
        fi

# make it so that we can learn the data format
# if we can find the word bam in the column then we do the bam path
# if we don't then we do the fastq path.

	
	# uncompress the files, the head pipe is for testing a subset of our data
	if [[ $D2 = *".bam" ]]; 
		then
			echo data is storred in unaligned bam format, copying and converting to fastq &>> $LOGDIR/$START/run.log
			rsync -av --no-p --no-o --no-g buckleyrm@lewis4-dtn.rnet.missouri.edu:$D1/$D2 $FQDIR/ &>> $LOGDIR/$START/run.log
			chmod uga+w $FQDIR/$D2

			# check if file transfered
			if [[ -e $FQDIR/$D2 ]]; 
				then
					echo file trasnfered &>> $LOGDIR/$START/run.log
				else
					echo transfer failed, exiting &>> $LOGDIR/$START/run.log
					exit
				fi

			samtools fastq --threads $THREADS -1 $FQDIR/$R1 -2 $FQDIR/$R2 $FQDIR/$D2 &>> $LOGDIR/$START/run.log

			# check if files successfully uncompressed
			if [[ -s $FQDIR/$R1 && -s $FQDIR/$R2 ]]; 
				then
					echo uncompressed read pair files found &>> $LOGDIR/$START/run.log
					rm $FQDIR/$D2
				else
					echo one or both uncompressed files are not found or are empty, exiting &>> $LOGDIR/$START/run.log
					exit
				fi

		else
			echo data is storred as compressed fastq, copying and uncompressing &>> $LOGDIR/$START/run.log
			rsync -av --no-p --no-o --no-g buckleyrm@lewis4-dtn.rnet.missouri.edu:$D1/$R1.gz $FQDIR/ &>> $LOGDIR/$START/run.log
			rsync -av --no-p --no-o --no-g buckleyrm@lewis4-dtn.rnet.missouri.edu:$D2/$R2.gz $FQDIR/ &>> $LOGDIR/$START/run.log
	
			chmod uga+w $FQDIR/$R1.gz
			chmod uga+w $FQDIR/$R2.gz

			# check if files trasnsfered
			if [[ -e $FQDIR/$R1.gz && -e $FQDIR/$R2.gz ]]; 
				then
					echo file trasnfered &>> $LOGDIR/$START/run.log
				else
					echo transfer failed, exiting &>> $LOGDIR/$START/run.log
					exit
				fi
			
			unpigz -p $THREADS $FQDIR/$R1.gz $FQDIR/$R2.gz
		fi


	echo begin mapping &>> $LOGDIR/$START/run.log
	# perform mapping
	(bwa mem -M -R $RG -t $THREADS $IDXDIR/$IDX $FQDIR/$R1 $FQDIR/$R2 | samtools view -Sb - > $MAPDIR/$SM.$LB.$IDfl.$IDln.bam) 2> $LOGDIR/$START/$SM.$LB.$IDfl.$IDln.aln.log

	#check for bam files and remove fastq files once reads are mapped
	if [ -s $MAPDIR/$SM.$LB.$IDfl.$IDln.bam ]
	then 
		echo bam file found, removing $FQDIR/{$R1,$R2} &>> $LOGDIR/$START/run.log
		rm -r $FQDIR/{$R1,$R2}
	else
		echo bam file not found or is empty, exiting &>> $LOGDIR/$START/run.log
		exit
	fi

	echo end mapping, begin sort &>> $LOGDIR/$START/run.log

	# sort and mark optical duplicates
	samtools sort --threads $THREADS -o $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam $MAPDIR/$SM.$LB.$IDfl.$IDln.bam &>> $LOGDIR/$START/run.log
	#gatk SortSam -I:$MAPDIR/$LB.$IDfl.$IDln.bam -O:$MAPDIR/$LB.$IDfl.$IDln.sorted.bam -SO:coordinate

	if [ -s $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam ]
	then
		echo sorted bam found, removing $MAPDIR/$SM.$LB.$IDfl.$IDln.bam &>> $LOGDIR/$START/run.log
		rm $MAPDIR/$SM.$LB.$IDfl.$IDln.bam
	else
		echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/run.log
		exit
	fi

	echo end sort, begin duplicate marking &>> $LOGDIR/$START/run.log

	gatk MarkDuplicates -I=$MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam -O=$MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.markedDup.bam -M=$MEDIR/$SM.$LB.$IDfl.$IDln.metrics 2> $LOGDIR/$START/$SM.$LB.$IDfl.$IDln.markDup.log
	# check for sorted bams and rm unsorted bams

	if [ -s $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.markedDup.bam ]
        then
                echo duplicate marked bam found, removing $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam &>> $LOGDIR/$START/run.log
                rm $MAPDIR/$SM.$LB.$IDfl.$IDln.sorted.bam
        else
                echo duplicate marked bam not found or is empty, exiting &>> $LOGDIR/$START/run.log
                exit
        fi
    echo end duplicate marking &>> $LOGDIR/$START/run.log
done


echo -e "read groups have all been processed, begin individual sample processing\n" &>> $LOGDIR/$START/run.log
# here we merge files, this is done on sample ID, which is the first name ID

for SAMPLE in $(awk '{print $5}' $SAMPLE_SHEET | uniq)
do
	echo -e "\n" begin merge for sample $SAMPLE "\n" &>> $LOGDIR/$START/run.log
	samtools merge --threads $THREADS $MAPDIR/$SAMPLE.bam $MAPDIR/$SAMPLE.*.markedDup.bam &>> $LOGDIR/$START/run.log
	if [ -s $MAPDIR/$SAMPLE.bam ]
	then 
		echo merged bam found, removing non-merged bams &>> $LOGDIR/$START/run.log
		rm $MAPDIR/$SAMPLE.*.markedDup.bam
	else
		echo merged bam not found or is empty, exiting &>> $LOGDIR/$START/run.log
		exit
	fi
	

	echo -e "\n" end merge for sample $SAMPLE, begin sort for sample $SAMPLE "\n" &>> $LOGDIR/$START/run.log
	# sort and mark optical duplicates
        samtools sort --threads $THREADS -o $MAPDIR/$SAMPLE.sorted.bam $MAPDIR/$SAMPLE.bam &>> $LOGDIR/$START/run.log
        #gatk SortSam -I:$MAPDIR/$LB.$IDfl.$IDln.bam -O:$MAPDIR/$LB.$IDfl.$IDln.sorted.bam -SO:coordinate

        if [ -s $MAPDIR/$SAMPLE.sorted.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SAMPLE.bam &>> $LOGDIR/$START/run.log
                rm $MAPDIR/$SAMPLE.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/run.log
                exit
        fi
    echo -e "\n" end sort for sample $SAMPLE, begin duplicate marking for sample $SAMPLE "\n" &>> $LOGDIR/$START/run.log
        gatk MarkDuplicates -I=$MAPDIR/$SAMPLE.sorted.bam -O=$MAPDIR/$SAMPLE.sorted.markedDup.bam -M=$MEDIR/$SAMPLE.metrics 2> $LOGDIR/$START/$SAMPLE.markDup.log
        # check for sorted bams and rm unsorted bams

        if [ -s $MAPDIR/$SAMPLE.sorted.markedDup.bam ]
        then
                echo sorted bam found, removing $MAPDIR/$SAMPLE.sorted.bam &>> $LOGDIR/$START/run.log
                rm $MAPDIR/$SAMPLE.sorted.bam
        else
                echo sorted bam not found or is empty, exiting &>> $LOGDIR/$START/run.log
                exit
        fi
    echo -e "\n" end duplicate marking for sample $SAMPLE "\n" &>> $LOGDIR/$START/run.log
done

echo mapping completed for sample sheet: $SAMPLE_SHEET &>> $LOGDIR/$START/run.log
cat $SAMPLE_SHEET &>> $LOGDIR/$START/run.log
echo done
