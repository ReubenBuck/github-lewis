#!/bin/bash
#SBATCH --partition BioCompute,Lewis  # for jobs < 2hrs try 'General'
#SBATCH -N1
#SBATCH -n51 # cores 
#SBATCH --mem 10G  # memory 
#SBATCH -t 2-00:00  # days-hours:minutes
#SBATCH --account=biocommunity  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=NA8EXQEX
#SBATCH --output=NA8EXQEX-%j.out  # %j is the unique jobID
#
## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send


# must run this from where the target files are located

OUT=NA8EXQEX_fullPath.result
OUT_DIR=/home/buckleyrm/storage.lyonslab/raw_data_backup
MD5_CHECK=/home/buckleyrm/storage.lyonslab/raw_data_backup/NA8EXQEX_fullPath.md5sum
THREADS=50
TMP_DIR=/home/buckleyrm/storage.lyonslab/TMP

touch $OUT_DIR/$OUT

N=$(expr $THREADS - 1)

for i in $(seq 1 $(wc -l $MD5_CHECK | cut -f1 -d" ")); do
	sleep 5s
	(	
		LINE=$(sed "${i}q;d" $MD5_CHECK)
		TMP_FILE=$(echo $LINE | cut -f1 -d" ")
		echo $LINE > $TMP_DIR/$TMP_FILE
		md5sum -c $TMP_DIR/$TMP_FILE >> $OUT_DIR/$OUT
		rm $TMP_DIR/$TMP_FILE
)&

if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
                wait
        fi

done

wait
