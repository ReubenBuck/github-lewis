#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition Lewis,BioCompute  # for jobs < 2hrs try 'General'
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MPI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=1 # cores per task
#SBATCH --mem-per-cpu=64G  # memory per core (default is 1GB/core)
#SBATCH --time 2-00:00  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=biocommunity  # investors will replace this with their account name
#
## labels and outputs
#SBATCH --job-name=sra_dl
#SBATCH --output=sra_dl-%A_%a.out  # %j is the unique jobID
#
## notifications
#SBATCH --mail-user=buckleyrm@missouri.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL  # which type of notifications to send
#
## array options
#
#-------------------------------------------------------------------------------

### NOTES



# invoke with RUN_SHEET=<sra run sheet with header> sbatch --array=2-$(wc -l < <sra run sheet with header>) dl_sra.sh

# set tmp dir for sra caching before starting, otherwise home will fill up and kill task


module load sratoolkit/sratoolkit-2.8.1-2


sleep $((RANDOM %10))

TMP=$(pwd)

if [ -d $TMP/.tmp_ncbi ]
then
        echo "tmp_ncbi already exists, continuing"
else
        echo "tmp_ncbi does not exist, creating"
        mkdir -p $TMP/.tmp_ncbi
        echo "/repository/user/main/public/root = \"$TMP/.tmp_ncbi\"" > $HOME/.ncbi/user-settings.mkfg
fi

#-------------------------------------------------------------------------------

echo "### Starting at: $(date) ###"

# random sleep time so runs don't collide
sleep $((RANDOM % 5))

# read sample name
ROW=$SLURM_ARRAY_TASK_ID
SM=$(cut -f 11 $RUN_SHEET | sed "${ROW}q;d")

# library name
LIB=$(cut -f 4 $RUN_SHEET | sed "${ROW}q;d")

# get run name
RUN=$(cut -f 9 $RUN_SHEET | sed "${ROW}q;d")

# is library paried or single
LAYOUT=$(cut -f 27 $RUN_SHEET | sed "${ROW}q;d")

# speices
SPECIES=$(cut -f 30 $RUN_SHEET | sed "${ROW}q;d")
SPECIES=${SPECIES// /_}


# create dir unless it already exists

if [ -d $SPECIES/$SM/$LIB ]
then
        echo "$SPECIES/$SM/$LIB dir exists, continuing"
else
        echo "$SPECIES/$SM/$LIB dir does not exist, creating"
        mkdir -p $SPECIES/$SM/$LIB
fi


echo -e "\nbegin fastq dump on sample $SM, library $LIB, srr $RUN \n"


# pull fastq RUN and store in ./sm/lib/
if [ $LAYOUT == 'PAIRED' ]
then
	fastq-dump --split-files --origfmt --gzip --outdir ./$SPECIES/$SM/$LIB $RUN
elif [ $LAYOUT == 'SINGLE' ]
then
	fastq-dump --origfmt --gzip --outdir ./$SPECIES/$SM/$LIB $RUN
else
	echo "check layout"
	exit
fi


#echo -e "\nfastq dump complete, renaming srr from $RUN to $LIB"
#mv ./$SM/$LIB/${RUN}_1.fastq.gz ./$SM/$LIB/${LIB}_1.fastq.gz
#mv ./$SM/$LIB/${RUN}_2.fastq.gz ./$SM/$LIB/${LIB}_2.fastq.gz

#echo -e "\nunzip read 1"
#gunzip ./$SM/$LIB/${LIB}_1.fastq.gz
#echo unzip done, split into runs
#~/scripts/bin/split_library_by_run ./$SM/$LIB/${LIB}_1.fastq
#echo split done, rm original
#rm ./$SM/$LIB/${LIB}_1.fastq
#echo zip new run fastq files
#gzip ./$SM/$LIB/${LIB}*_1.fastq

#echo -e "\nunzip read 2"
#gunzip ./$SM/$LIB/${LIB}_2.fastq.gz
#echo unzip done, split into runs
#~/scripts/bin/split_library_by_run ./$SM/$LIB/${LIB}_2.fastq
#echo split done, rm original
#rm ./$SM/$LIB/${LIB}_2.fastq
#echo zip new fun fastq files
#gzip ./$SM/$LIB/${LIB}*_2.fastq


echo "### Ending at: $(date) ###" #!/bin/bash
