Using the shotgun-sunbeam pipeline for shotgun metagenomics sequencing
======================================================================

TBD

## Working on the cluster

Here we provide basic commands for connecting and submitting jobs on the cluster. See https://docs.massive.org.au/ for more information.

```
# On Mac
ssh [username]@m3.massive.org.au

# On Linux
ssh -l [username] m3.massive.org.au

# Start a new smux session (n) for example with 20 cpus (--ntasks) for 1 day (--time)
smux n --ntasks=20 --time=1-00:00:00

# List ongoing smux sessions
smux l

# List ongoing jobs
show_job

# Load modules
module load [module]

# Cancel a smux session
scancel [JOBID]
```

## Installing Sunbeam on the cluster (to do only ONCE)

Sunbeam is a pipeline written in snakemake that simplifies and automates many of the steps in metagenomic sequencing analysis. See https://sunbeam.readthedocs.io/en/latest/index.html for more information.

```
module load git
git clone -b stable https://github.com/sunbeam-labs/sunbeam sunbeam-stable
cd sunbeam-stable
chmod u+w .git/objects/pack/* # Very important, to work around permission issues on the cluster!
bash install.sh
```

You will also need some extensions.

```
module load git
source activate sunbeam
cd sunbeam-stable

# Kraken 2 extension
git clone https://github.com/louiejtaylor/sbx_kraken2/ extensions/sbx_kraken2
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

1) Blast databases for nucleic acid (nt) and protein (nr) mapping module load blast

```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
gunzip nt.gz
makeblastdb -in nt -out nt -dbtype nucl

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
makeblastdb -in nr -out nr -dbtype prot
```

2) Kraken databases for taxonomy

```
# Kraken
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
tar -xvzf minikraken_20171019_8GB.tgz

# Kraken2
wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
tar -xvzf minikraken_8GB_202003.tgz
Another database (eg fungi)
Check https://github.com/mw55309/Kraken_db_install_scripts for more details of construction a Kraken-compatible database for fungi or other. Metaphlan extension contains an eucaryotic database as well.
```

3) Other databases 

Check https://github.com/mw55309/Kraken_db_install_scripts for more details of construction a Kraken-compatible database for fungi or other. Metaphlan extension contains an eucaryotic database as well.

## Initialise your sunbeam project

sunbeam init takes one required argument: a path to your project folder. This folder will be created if it doesn’t exist. You can also specify the path to your gzipped fastq files, and Sunbeam will try to guess how your samples are named, and whether they’re paired.

```
sunbeam init --data_fp [/path/to/fastq/files /path/to/my_project]
```

Because the default sunbeam_config.yml does not contain the extensions parameters, update it by running:

```
cat extensions/sbx_kraken2/config.yml >> /path/to/my_project/sunbeam_config.yml
cat extensions/sbx_subsample/config.yml >> /path/to/my_project/sunbeam_config.yml
cat extensions/sbx_metaphlan/config.yml >> /path/to/my_project/sunbeam_config.yml
cat extensions/sbx_eggnog/config.yml >> /path/to/my_project/sunbeam_config.yml
```

In your project directory directory, a new config file and a new sample list were created (by default named sunbeam_config.yml and samplelist.csv, respectively). Edit the config file in your favorite text editor and samplelist.csv if necessary. You may want to check the paths to your project, databases, adapter sequences etc. An example of the sunbeam_config.yml is provided (here)[git].

## Citation

If you used this repository in a publication, please mention its url.

In addition, you may cite the tools used by this pipeline:

* **Sunbeam:** EL Clarke, LJ Taylor, C Zhao et al. Sunbeam: an extensible pipeline for analyzing metagenomic sequencing experiments. Microbiome 7:46 (2019).

## Rights

* Copyright (c) 2020 Respiratory Immunology lab, Monash University, Melbourne, Australia.
* License: The R Notebook template (.Rmd) is provided under the MIT license (See LICENSE.txt for details)
* Authors: C. Pattaroni, B.J. Marsland
