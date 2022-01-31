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

### Installing extensions

You will also need some extensions (in the `extensions` folder of sunbeam).

```
# Kraken2bracken extension
git clone https://github.com/respiratory-immunology-lab/sbx_kraken2bracken sunbeam-stable/extensions/sbx_kraken2bracken
```

Because the default sunbeam_config.yml does not contain the extensions parameters, update it by running:

```
cat sunbeam-stable/extensions/sbx_kraken2/config.yml >> /path/to/my_project/sunbeam_config.yml
```

## Databases

1) Host genome(s) for decontamination
2) Kraken databases for taxonomy
3) Bracken databases (related to kraken2)

All these databases are *available on the cluster* at `~/of33/Databases/shotgun`. Don't re-build your own unless absolutely necessary. If necessary, follow the instructions provided [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/tree/master/databases).

## Initialise your sunbeam project

sunbeam init takes one required argument: a path to your project folder. This folder will be created if it doesn’t exist. You can also specify the path to your gzipped fastq files, and Sunbeam will try to guess how your samples are named, and whether they’re paired.

```
source activate sunbeam
sunbeam init --data_fp /path/to/fastq/files /path/to/my_project
```

In your project directory directory, a new config file and a new sample list were created (by default named sunbeam_config.yml and samplelist.csv, respectively). Edit the config file in your favorite text editor and samplelist.csv if necessary. You may want to check the paths to your project, databases, adapter sequences etc. An example of the sunbeam_config.yml is provided [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/Sunbeam/blob/master/sunbeam_config.yml).

## Running sunbeam

These are the parameters for using the `--partition=genomics --qos=genomics` partition on the cluster.

```
# Filtering and host sequences decontamination
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=sunbeam_all_decontam --account=of33 --time=04:00:00 --mem-per-cpu=8G --ntasks=1 --cpus-per-task=16 --partition=genomics --qos=genomics" -j 30 -w 60 -p all_decontam

# Classification with Kraken2 and Bracken correction
sunbeam run --configfile sunbeam_config.yml --cluster "sbatch --job-name=sunbeam_all_kraken2bracken --account=of33 --time=04:00:00 --mem-per-cpu=8G --ntasks=1 --cpus-per-task=16 --partition=genomics --qos=genomics" -j 30 -w 60 -p --use-conda all_kraken2bracken
```

## Getting host versus non-host read counts

The total number of reads and the non-mapping number of reads can be retrieved using this wrapper script.

```
for f in sunbeam_output/qc/cleaned/*_1.fastq.gz;
do basename=${f%_1*};
sample=${basename#*/*/*/*};
clean=sunbeam_output/qc/cleaned/${sample}_1.fastq.gz;
decontam=sunbeam_output/qc/decontam/${sample}_1.fastq.gz;
cleancount=`gunzip -c $clean | awk 'END{print(NR/4)}'`
decontamcount=`gunzip -c $decontam | awk 'END{print(NR/4)}'`
echo $sample$'\t'$cleancount$'\t'$decontamcount;
done > sunbeam_output/qc/counts.tsv
```

## Functional profiling with MetaLAFFA

Functional profiling is performed using the MetaLAFFA tool (more information available [here](https://github.com/borenstein-lab/MetaLAFFA)). This tool is also written in Snakemake but a few steps are required to enable it to run on the cluster.

### Installation

As MetaLAFFA requires a specific version of python, the conda environment needs to be as following:

```
conda create -n metalaffa python=3.6.10 
conda activate metalaffa
conda install -c bioconda -c borenstein-lab metalaffa
```

If you see any issue with the installation, try using `conda clean --all` and `conda config --remove channels conda-forge` commands before you create the new conda environment.

### Databases installation

MetaLAFFA automatically installs the databases in the conda environment path. To choose another installation path, change the `file_organisation.py` file located in `~/miniconda3/envs/metalaffa/lib/python3.6/config/` with a new database location. 

Note that all databases are already installed in the cluster's database folder `~/of33/Databases/shotgun/metalaffa/`. See an example of the updated `file_organisation.py` file [here](https://github.com/respiratory-immunology-lab/microbiome-shotgun/blob/master/MetaLAFFA/file_organisation.py).

If you still need to install the databases, run:

```
# Download and prepare default reference databases
prepare_databases.py -hr -km -u -c
```

### Creation a new metalaffa project folder

With your MetaLAFFA conda environment active, you can create a new metalaffa project directory as following. This will also copy some of the parameter files needed to run the pipeline.

```
# Create a metalaffa project directory
create_new_MetaLAFFA_project.py [location_of_your_choice]/metalaffa

# Enter the directory
cd metalaffa
```

### Change configuration files to run MetaLAFFA on the cluster

A few files need to be updated or added to enable the use of MetaLAFFA on the cluster and tell MetaLAFFA to skip the QC and host decontamination steps (already performed by sunbeam).

1) `cluster.py` located in `[location_of_your_choice]/metalaffa/config/` needs to be updated to [THIS](https://github.com/respiratory-immunology-lab/microbiome-shotgun/blob/master/MetaLAFFA/cluster.py).
2) `pipeline_steps.txt` located in `[location_of_your_choice]/metalaffa/` needs to be updated to [THIS](https://github.com/respiratory-immunology-lab/microbiome-shotgun/blob/master/MetaLAFFA/pipeline_steps.txt).
3) A new job submission file specific to M3 needs `m3_submission_wrapper.py` needs to be added in `[location_of_your_choice]/src/`, this file is available [HERE](https://github.com/respiratory-immunology-lab/microbiome-shotgun/blob/master/MetaLAFFA/m3_submission_wrapper.py).

Add executing permission to the `m3_submission_wrapper.py` file:

```
chmod +x [location_of_your_choice]/metalaffa/src/m3_submission_wrapper.py
```

### Running MetaLAFFA

Copy your clean and host decontaminated files into the MetaLAFFA `data` folder.

```
cp [sunbeam_location]/sunbeam_output/qc/decontam/*fastq.gz [location_of_your_choice]/metalaffa/data/
```

MetaLAFFA requires specific sample names finishing in `.R1.fastq.gz` and `.R2.fastq.gz`. If required, rename files as following:

```
# If files finish with _1.fastq.gz and _2.fastq.gz
rename _1. .R1. *
rename _2. .R2. *

# If files finish with _R1.fastq.gz and _R2.fastq.gz
rename _R1. .R1. *
rename _R2. .R2. *
```

You are now finally ready to use MetaLAFFA. Don't forget the `--use-cluster` flag to run the pipeline on the cluster.

```
# Run MetaLAFFA (in your metalaffa directory)
./MetaLAFFA.py --use-cluster
```

## Citation

If you used this repository in a publication, please mention its url.

In addition, you may cite the tools used by this pipeline:

* **Sunbeam:** EL Clarke, LJ Taylor, C Zhao et al. Sunbeam: an extensible pipeline for analyzing metagenomic sequencing experiments. Microbiome 7:46 (2019).
* **MetaLAFFA:** Eng, A., Verster, A. J., & Borenstein, E. (2020). MetaLAFFA: a flexible, end-to-end, distributed computing-compatible metagenomic functional annotation pipeline. BMC bioinformatics, 21(1), 1-9.

## Rights

* Copyright (c) 2021 Respiratory Immunology lab, Monash University, Melbourne, Australia.
* License: The R Notebook template (.Rmd) is provided under the MIT license (See LICENSE.txt for details)
* Authors: C. Pattaroni
