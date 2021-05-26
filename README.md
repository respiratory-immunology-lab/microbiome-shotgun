Using the shotgun-sunbeam pipeline for shotgun metagenomics sequencing
======================================================================

This pipeline is based on Sunbeam pipeline, which provides a foundation to build more in-depth analyses and to enable comparisons in metagenomic sequencing experiments by removing problematic, low-complexity reads and standardizing post-processing and analytical steps. Sunbeam is written in Python using the Snakemake workflow management software. This pipeline has been optimised to work on M3 cluster.

## Working on the cluster

Click [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/tree/master/cluster) for more information about how to work on the cluster and download data directly from basespace.

## Installing Sunbeam on the cluster

Sunbeam is a pipeline written in snakemake that simplifies and automates many of the steps in metagenomic sequencing analysis. See https://sunbeam.readthedocs.io/en/latest/index.html for more information.

```
# Load git module on cluster
module load git

# Clone sunbeam
git clone -b stable https://github.com/sunbeam-labs/sunbeam sunbeam-stable

# Change permissions (if not, creates issues with cluster)
cd sunbeam-stable
chmod u+w .git/objects/pack/*

# Install sunbeam
bash install.sh

# Test the installation
tests/run_tests.bash -e sunbeam
```

You will also need some extensions.

```
# Kraken2 extension
git clone https://github.com/cpattaroni/sbx_kraken2 extensions/sbx_kraken2
cat extensions/sbx_kraken2/config.yml >> /path/to/my_project/sunbeam_config.yml

# Braken extension
git clone https://github.com/cpattaroni/sbx_bracken extensions/sbx_braken
cat extensions/sbx_braken/config.yml >> /path/to/my_project/sunbeam_config.yml
```

## Databases

1) Host genome(s) for decontamination
2) Kraken databases for taxonomy
3) Bracken databases (related to kraken2)
4) Blast databases for nucleic acid (nt) and protein (nr) mapping

All these databases are *available on the cluster* at `~/of33/Databases/shotgun`. Don't re-build your own unless absolutely necessary. If necessary, follow the instructions provided [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/tree/master/databases).

## Initialise your sunbeam project

sunbeam init takes one required argument: a path to your project folder. This folder will be created if it doesn’t exist. You can also specify the path to your gzipped fastq files, and Sunbeam will try to guess how your samples are named, and whether they’re paired.

```
source activate sunbeam
sunbeam init --data_fp [/path/to/fastq/files /path/to/my_project]
```

Because the default sunbeam_config.yml does not contain the extensions parameters, update it by running:

```
cat extensions/sbx_kraken2/config.yml >> /path/to/my_project/sunbeam_config.yml
```

In your project directory directory, a new config file and a new sample list were created (by default named sunbeam_config.yml and samplelist.csv, respectively). Edit the config file in your favorite text editor and samplelist.csv if necessary. You may want to check the paths to your project, databases, adapter sequences etc. An example of the sunbeam_config.yml is provided [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/blob/master/sunbeam_config.yml).

## Running sunbeam

These are the parameters for using the `--partition=genomics --qos=genomics` partition on the cluster.

```
# Filtering and host sequences decontamination
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=all_decontam --account=of33 --time=04:00:00 --mem=50G --ntasks=24 --partition=genomics --qos=genomics" -j 48 -w 30 -p all_decontam

# Classification with Kraken2
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=all_kraken2 --account=of33 --time=04:00:00 --mem=50G --ntasks=24 --partition=genomics --qos=genomics" -j 48 -w 30 -p --use-conda all_kraken2

# Assembly with Megahit
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=all_assembly --account=of33 --time=04:00:00 --mem=50G --ntasks=24 --partition=genomics --qos=genomics" -j 48 -w 30 -p all_assembly

# Annotation
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=all_annotate --account=of33 --time=04:00:00 --mem=50G --ntasks=24 --partition=genomics --qos=genomics" -j 48 -w 30 -p all_annotate

# Functional annotation
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=all_eggnog_bac --account=of33 --time=04:00:00 --mem=50G --ntasks=24 --partition=genomics --qos=genomics" -j 48 -w 30 -p --use-conda all_eggnog_bac
```

## Citation

If you used this repository in a publication, please mention its url.

In addition, you may cite the tools used by this pipeline:

* **Sunbeam:** EL Clarke, LJ Taylor, C Zhao et al. Sunbeam: an extensible pipeline for analyzing metagenomic sequencing experiments. Microbiome 7:46 (2019).

## Rights

* Copyright (c) 2021 Respiratory Immunology lab, Monash University, Melbourne, Australia.
* License: The R Notebook template (.Rmd) is provided under the MIT license (See LICENSE.txt for details)
* Authors: C. Pattaroni
