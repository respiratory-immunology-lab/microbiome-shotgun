#!/bin/bash

# Set the location of the MEGAHIT folder
MEGAHIT_FP="./sunbeam_output/assembly/megahit"

rm -r megahit_contig_stats.csv # Remove the existing file if it exists
touch megahit_contig_stats.csv # Replace with a new file

# Make a header
echo 'sample,contigs,total,min,max,avg,N50' >> megahit_contig_stats.csv

for sample in `ls $MEGAHIT_FP | awk -F'[_]' '{print $1 "_" $2}'` # Loops through the sample names
# This works for samples named like 'patientID_S1' or 'uniqueID_SampleNumber' etc.
# E.g. '180s2_S32'
# Change this if necessary
do
  # Prepare a variable for the filename - subfolders in the megahit folder have the '_asm' suffix
	file=$MEGAHIT_FP/$sample*_asm*/log
  # Print comma-separated values and add to the file
	cat $file | grep STAT | awk -v sample=$sample '{print sample "," $3 "," $6 "," $9 "," $12 "," $15 "," $18}' >> \
	megahit_contig_stats.csv
done
