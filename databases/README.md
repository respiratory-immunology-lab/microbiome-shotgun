# Databases building


## Recommended smux parameters for database building

```
smux n --nodes=1 --ntasks=1 --cpuspertask=36 --mem=100000 -J Databases --time=5-00:00:00
```

## Host genome(s)

We need the host genomes to remove host reads before metagenomics analysis. Of note, these need to be located in a separate folder, be decompressed, and end up with `.fasta`. 

```
# Mouse
wget http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz
mv Mus_musculus.GRCm38.dna_sm.primary_assembly.fa Mus_musculus.GRCm38.dna_sm.primary_assembly.fasta

# Human
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa Homo_sapiens.GRCh38.dna_sm.primary_assembly.fasta
```

## Blast databases for nucleic acid (nt) and protein (nr) mapping module load blast

```
# Load blast
module load blast

# Nt database
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip -d nt.gz
makeblastdb -in nt -out nt -dbtype nucl

# Nr database
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip -d nr.gz
makeblastdb -in nr -out nr -dbtype prot
```

## Blast database of virulence factors

```
# Load blast
module load blast

# Nt database
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
gunzip -d VFDB_setB_nt.fas.gz
makeblastdb -in VFDB_setB_nt.fas -out VFDB_setB_nt -dbtype nucl

# Nr database
wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
gunzip -d VFDB_setB_pro.fas.gz

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
kraken2-build --download-library archaea --db [mydatabase]
kraken2-build --download-library bacteria --db [mydatabase]
kraken2-build --download-library fungi --db [mydatabase]
kraken2-build --download-library viral --db [mydatabase]

# Build the database (takes time)
kraken2-build --build --db [mydatabase]
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

