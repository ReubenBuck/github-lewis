#!/bin/bash

module load samtools

module load bcftools

EXCDIR=/home/buckleyrm/software/excord
GIGDIR=/home/buckleyrm/software/giggle/bin
STIDIR=/home/buckleyrm/software/stix/bin

BAMDIR=/home/buckleyrm/storage.lyonslab/results/Felis_catus/bams

REF=/storage/htc/lyonslab/cat_ref/Felis_catus_9.0.fa

OUTDIRDIS=/storage/hpc/group/UMAG/WORKING/buckleyrm/dwarfism/all/discord

OUTDIRIDX=/storage/hpc/group/UMAG/WORKING/buckleyrm/dwarfism/all/stix_idx

PED=/home/buckleyrm/dwf.all.ped

REGION="chrB1:174874754-174892181"

# create out dirs
if [[ ! -d $OUTDIRDIS ]]; then
	mkdir -p $OUTDIRDIS
fi

if [[ ! -d $OUTDIRIDX ]]; then
	mkdir -p $OUTDIRIDX
fi

# extract discordant reads

ROWS=$(seq 1 $(tail -n+2 $PED | wc -l))

BAMS=($(cut -f2 $PED))
BEDS=($(cut -f3 $PED))

for i in $ROWS; do
	#sleep $((RANDOM % 2))
	if [[ ! -s $OUTDIRDIS/${BEDS[$i]} ]]; then
		(
			samtools view -b $BAMDIR/${BAMS[$i]/.sort.markDup.realign.bam/}/${BAMS[$i]} $REGION \
			| $EXCDIR/excord \
			--discordantdistance 500 \
			--fasta $REF \
			/dev/stdin \
			| LC_ALL=C sort --buffer-size 2G -k1,1 -k2,2n -k3,3n \
			| bgzip -c > $OUTDIRDIS/${BEDS[$i]}
			) 
	fi

done

wait

# index
if [[ ! -z $OUTDIRIDX ]]; then
	$GIGDIR/giggle index -i "$OUTDIRDIS/*.bed.gz" -o $OUTDIRIDX -s -f
	# feeding beds needs work, not printing right	
fi

#create db
if [[ ! -s $PED.db ]]; then
	$STIDIR/stix -i $OUTDIRIDX -p $PED -c 3 -d $PED.db
fi
