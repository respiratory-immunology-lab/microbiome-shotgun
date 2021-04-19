#!/bin/bash

# Create a new directory for the fastq files
mkdir rawfastq

# Download data from basespace
./bs -c Australia download project -i [projectnumber] -o basespace_output --extension=fastq.gz

# Assuming we have the BaseSpace project downloaded in rawfastq directory
find basespace_output -name "*.fastq.gz" -exec cp {} rawfastq \;

# Remove initial folder to make space
rm -rf basespace_output

# Merge lanes 1 and 2
cd rawfastq
for f in *.fastq.gz
  do
  Basename=${f%_L00*}
  ## merge R1
  ls ${Basename}_L00*_R1_001.fastq.gz | xargs cat > ${Basename}_merged_R1.fastq.gz
  ## merge R2
  ls ${Basename}_L00*_R2_001.fastq.gz | xargs cat > ${Basename}_merged_R2.fastq.gz
  done

# Check whether the merging has worked
#-bash-4.2$ zcat []_L001_R1_001.fastq.gz | echo $((`wc -l`/4))
#-bash-4.2$ zcat []_L002_R1_001.fastq.gz | echo $((`wc -l`/4))
#-bash-4.2$ zcat []_merged_R1.fastq.gz | echo $((`wc -l`/4))

# Remove individual files to make space
rm -rf *L00*

# Save files to Google Drive (make sure you it's configured for your account)
module load rclone
rclone copy rawfastq gdrive:/02_Data/Australia/NovaSeq

