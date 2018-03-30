#!/bin/bash

for file in $(ls *.bam); do
echo $(cat "$file".md5)	$file
done
