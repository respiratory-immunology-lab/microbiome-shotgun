# Basespace

Basespace is the Illumina tool enabling the download of sequencing data.

```
# Create a new directory for the fastq files
mkdir rawfastq

# Download data from basespace
./bs -c Australia download project -i [projectnumber] -o basespace_output --extension=fastq.gz

# Assuming we have the BaseSpace project downloaded in rawfastq directory
find basespace_output -name "*.fastq.gz" -exec cp {} rawfastq \;

# Remove initial folder to make space
rm -rf basespace_output
```

# Lanes merging

Sometimes, libraries are sequenced on different lanes. The lines below can be used to merge files from different lanes.

```
# Merge lanes 1 and 2
cd rawfastq
for f in *.fastq.gz
  do
  Basename=${f%_L00*}
  ## merge R1
  ls ${Basename}_L00*_R1_001.fastq.gz | xargs cat > ${Basename}_R1.fastq.gz
  ## merge R2
  ls ${Basename}_L00*_R2_001.fastq.gz | xargs cat > ${Basename}_R2.fastq.gz
  done

# Remove individual files to make space
rm -rf *L00*
```

# Use of rclone

The data can be stored to Google Drive using rclone tool.

```
# Save files to Google Drive (make sure you it's configured for your account)
module load rclone
rclone copy -v --fast-list --max-backlog=999999 --drive-chunk-size=512M rawfastq gdrive:/02_Data/MonashData/NovaSeq/NovaSeqXX/Data
```
