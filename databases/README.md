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

## Blast databases for nucleic acid (nt) and protein (nr) mapping module load blast

```
# Activate sunbeam
source activate sunbeam

# Nt database
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip -d nt.gz
makeblastdb -in nt -out nt -dbtype nucl

# Nr database
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip -d nr.gz
makeblastdb -in nr -out nr -dbtype prot
```

## Blast database of virulence factors (optional)

```
# Activate sunbeam
source activate sunbeam

# Nt database
wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
gunzip -d VFDB_setB_nt.fas.gz
makeblastdb -in VFDB_setB_nt.fas -out VFDB_setB_nt -dbtype nucl

# Nr database
wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
gunzip -d VFDB_setB_pro.fas.gz
makeblastdb -in VFDB_setB_pro.fas -out VFDB_setB_pro -dbtype prot
```

## Kraken databases for taxonomy

There are some small pre-compiled databases available. However, because we want to go in more depth and are interested in fungi as well we will build our own.

```
# Install kraken2
git clone https://github.com/DerrickWood/kraken2 [mydirectory]
cd [mydirectory]
bash install_kraken2.sh

# Install the taxonomy
kraken2-build --download-taxonomy --db [mydatabase]

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

