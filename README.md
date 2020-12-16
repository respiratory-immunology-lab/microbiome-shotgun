Using the shotgun-sunbeam pipeline for shotgun metagenomics sequencing
======================================================================

This pipeline is based on Sunbeam pipeline, which provides a foundation to build more in-depth analyses and to enable comparisons in metagenomic sequencing experiments by removing problematic, low-complexity reads and standardizing post-processing and analytical steps. Sunbeam is written in Python using the Snakemake workflow management software. This pipeline has been optimised to work on M3 cluster.

## Working on the cluster

Here we provide basic commands for connecting and submitting jobs on the cluster. See https://docs.massive.org.au/ for more information on the cluster.

```
# Connect to cluster
ssh [username]@m3.massive.org.au # Mac users
ssh -l [username] m3.massive.org.au # Linux users

# Start a new smux session (n) for example with 20 cpus (--ntasks) for 1 day (--time)
smux n --ntasks=20 --time=1-00:00:00

# List ongoing smux sessions
smux l

# List ongoing jobs
show_job

# Load modules
module load [module]

# Cancel a smux session
scancel [jobID]
```

## Downloading ata from BaseSpace

See https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview for downloading data from BaseSpace straight to the cluster. 

```
# Download data from BaseSpace
./bs -c Australia download project -i [projectID] -o [directory]
```

## Installing Sunbeam on the cluster

Sunbeam is a pipeline written in snakemake that simplifies and automates many of the steps in metagenomic sequencing analysis. See https://sunbeam.readthedocs.io/en/latest/index.html for more information.

```
# Load git mocule on cluster
module load git

# Clone sunbeam
git clone -b stable https://github.com/sunbeam-labs/sunbeam sunbeam-stable

# Change permissions (if not, creates issues with cluster)
cd sunbeam-stable
chmod u+w .git/objects/pack/*

# Install sunbeam
bash install.sh

# Test the installation
cd sunbeam-stable
tests/run_tests.bash -e sunbeam
```

You will also need some extensions.

```
module load git
source activate sunbeam
cd sunbeam-stable

# Kraken 2 extension
git clone https://github.com/louiejtaylor/sbx_kraken2 extensions/sbx_kraken2
cat extensions/sbx_kraken2/config.yml >> /home/cpat0003/of33_scratch/Shotgun/MD4_project/sunbeam_config.yml

# Subsample extension (for assembly)
git clone https://github.com/sunbeam-labs/sbx_subsample/ extensions/sbx_subsample
cat extensions/sbx_subsample/config.yml >> /home/cpat0003/of33_scratch/Shotgun/MD4_project/sunbeam_config.yml

# eggNOG (functionnal annotation)
git clone https://github.com/ArwaAbbas/sbx_eggnog/ extensions/sbx_eggnog
cat extensions/sbx_eggnog/config.yml >> /home/cpat0003/of33_scratch/Shotgun/MD4_project/sunbeam_config.yml

# sbx_humann
...Coming soon...
```

## Databases

1) Host genome(s)

We need the host genomes to remove host reads before metagenomics analysis. Of note, these need to be located in a separate folder, be decompressed, and end up with `.fasta`. 

```
# Mouse
wget http://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculuss.GRCm38.dna_sm.primary_assembly.fa.gz

# Human
wget http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
```

2) Blast databases for nucleic acid (nt) and protein (nr) mapping module load blast

```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz
makeblastdb -in nt -out nt -dbtype nucl

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
makeblastdb -in nr -out nr -dbtype prot
```

3) Kraken databases for taxonomy

There are some small pre-compiled databases available. However, because we want to go in more depth and are interested in fungi as well we will build our own. A wrapper script for the database construction `builddatabase.sh` is provided [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/builddatabase.sh).

```
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

## Initialise your sunbeam project

sunbeam init takes one required argument: a path to your project folder. This folder will be created if it doesn’t exist. You can also specify the path to your gzipped fastq files, and Sunbeam will try to guess how your samples are named, and whether they’re paired.

```
source activate sunbeam
sunbeam init --data_fp [/path/to/fastq/files /path/to/my_project]
```

Because the default sunbeam_config.yml does not contain the extensions parameters, update it by running:

```
cat extensions/sbx_kraken2/config.yml >> /path/to/my_project/sunbeam_config.yml
cat extensions/sbx_subsample/config.yml >> /path/to/my_project/sunbeam_config.yml
cat extensions/sbx_metaphlan/config.yml >> /path/to/my_project/sunbeam_config.yml
cat extensions/sbx_eggnog/config.yml >> /path/to/my_project/sunbeam_config.yml
```

In your project directory directory, a new config file and a new sample list were created (by default named sunbeam_config.yml and samplelist.csv, respectively). Edit the config file in your favorite text editor and samplelist.csv if necessary. You may want to check the paths to your project, databases, adapter sequences etc. An example of the sunbeam_config.yml is provided [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/).

## Citation

If you used this repository in a publication, please mention its url.

In addition, you may cite the tools used by this pipeline:

* **Sunbeam:** EL Clarke, LJ Taylor, C Zhao et al. Sunbeam: an extensible pipeline for analyzing metagenomic sequencing experiments. Microbiome 7:46 (2019).

## Rights

* Copyright (c) 2020 Respiratory Immunology lab, Monash University, Melbourne, Australia.
* License: The R Notebook template (.Rmd) is provided under the MIT license (See LICENSE.txt for details)
* Authors: C. Pattaroni, B.J. Marsland
