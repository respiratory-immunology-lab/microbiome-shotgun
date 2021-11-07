# Databases building


## Recommended smux parameters for database building

```
smux n --nodes=1 --ntasks=1 --mem=10G -J Databases --time=7-00:00:00
```

## Host genome(s)

We need the host genomes to remove host reads before metagenomics analysis. Of note, these need to be located in a separate folder, be decompressed, and end up with `.fasta`. For human data, the best is to use a decoy version of the genome.

```
# Mouse
wget -c http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
gzip -d Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
mv Mus_musculus.GRCm38.dna_sm.primary_assembly.fa Mus_musculus.GRCm38.dna_sm.primary_assembly.fasta

# Human
wget -c http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_re                                                ference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
mv hs37d5.fa hs37d5.fasta
```

## Host and rRNA tRNA reference files (for metatranscriptomics)

rRNA and tRNA files can be found [here](https://github.com/elfrouin/transcriptM/tree/master/databases/1-SortMeRNA). 

```
# Human transcriptome files
wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gzip -d Homo_sapiens.GRCh38.cds.all.fa.gz
wget -c http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
gzip -d Homo_sapiens.GRCh38.ncrna.fa.gz

# Other eucaryotic and procaryotic files can be found in the link above
```

## Kraken databases for taxonomy

There are some small pre-compiled databases available. However, because we want to go in more depth and are interested in fungi as well we will build our own.

```
# Install kraken2
git clone https://github.com/DerrickWood/kraken2 [mydirectory]
cd [mydirectory]
bash install_kraken2.sh

# Install the taxonomy
~/kraken2/kraken2-build --download-taxonomy --db [mydatabase]

# Load blast module (for low complexity sequences masking)
module load blast

# Download reference libraries
~/kraken2/kraken2-build --download-library archaea --db [mydatabase]
~/kraken2/kraken2-build --download-library bacteria --db [mydatabase]
~/kraken2/kraken2-build --download-library fungi --db [mydatabase]
~/kraken2/kraken2-build --download-library viral --db [mydatabase]

# Build the database (takes time)
~/kraken2/kraken2-build --build --db [mydatabase]
```

## Braken for corrected species abundances

```
# Install braken
git clone https://github.com/jenniferlu717/Bracken [mydirectory]
cd [mydirectory]
bash install_bracken.sh

# Build the database
bracken-build -d [mydatabase] -t [numberofthreads] -x [directorywherekraken2isinstalled]
```

## RefSeq database for annotation (proteins)

```
# Activate sunbeam
source activate sunbeam

# Get the RefSeq summary files for viral, bacterial and fungal genomes
curl -o viral_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
curl -o bacteria_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
curl -o fungi_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt

# Create a single file with links for downloads (protein coding fasta files of complete genomes)
cat viral_summary.txt bacteria_summary.txt fungi_summary.txt > assembly_summary.txt
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths
awk 'BEGIN{FS=OFS="/";filesuffix="protein.faa.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths > ftpfilepaths 

# Download all genomes and create a single fasta file
wget -i ftpfilepaths #35009 genomes downloaded Nov 2021
gunzip *.gz
cat *.faa > refseq.fa
rm -rf *_protein.faa

# Create database
makeblastdb -in refseq.fa -title refseq -dbtype prot -hash_index -max_file_sz '4GB'
```

## CARD database for resistance genes (optional)

```
# Activate sunbeam
source activate sunbeam

# Create database
curl -o broadstreet-v3.1.4.tar.bz2 https://card.mcmaster.ca/download/0/broadstreet-v3.1.4.tar.bz2
tar -xf broadstreet-v3.1.4.tar.bz2
makeblastdb -in protein_fasta_protein_homolog_model.fasta -title card_protein -dbtype prot -hash_index -max_file_sz '4GB'
makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta -title card_nucl -dbtype nucl -hash_index -max_file_sz '4GB'
```

## Virulence factor database (optional)

```
# Activate sunbeam
source activate sunbeam

# Virulence factors database
wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
gunzip -d VFDB_setB_nt.fas.gz
makeblastdb -in VFDB_setB_nt.fas -title VFDBnucl -dbtype nucl -hash_index -max_file_sz '4GB'

# Virulence factors database
wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
gunzip -d VFDB_setB_pro.fas.gz
makeblastdb -in VFDB_setB_pro.fas -title VFDBprotein -dbtype prot -hash_index -max_file_sz '4GB'
```
