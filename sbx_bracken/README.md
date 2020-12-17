# sbx_bracken
This is an extension to the [Sunbeam pipeline](https://github.com/sunbeam-labs/sunbeam) to run Bracken (Bayesian Reestimation of Abundance with KrakEN). Braken uses the taxonomy labels assigned by Kraken, a highly accurate metagenomics classification algorithm, to estimate the number of reads originating from each species present in a sample. See [Braken](https://ccb.jhu.edu/software/bracken/index.shtml) for more information.

## Installing

Clone the repo into your sunbeam `extensions/` folder and add the new options to your existing configuration file. Make sure you've [installed Sunbeam](https://sunbeam.readthedocs.io/en/latest/quickstart.html) first!

```
source activate sunbeam
cd $SUNBEAM_DIR
git clone https://github.com/respiratory-immunology-lab/microbiome-shotgun/sbx_bracken extensions/sbx_bracken
```
Add the options to your config file (replace "sunbeam_config.yml" with the name of your config file).

```
cat extensions/sbx_kraken2/config.yml >> sunbeam_config.yml
```

## Configuration

In your configuration file, make sure you choose how many threads you'd like kraken2 to use, as well as providing the path to your database:

```
sbx_kraken2:
  threads: 4
  db_file: '/path/to/your/bracken/database_dir'
```

## Running

Finally, run Sunbeam as usual (with the `--use-conda` flag) with your extension's target rule specified:

```
sunbeam run --use-conda --configfile=sunbeam_config.yml all_bracken
```

The `--use-conda` flag is crucial as it enables running bracken in an isolated environment.

This rule generates output in the following location: `sunbeam_output/classify/bracken/`

## Contents

 - `sbx_bracken_env.yml` specifies the extension's dependencies and provides the environment that Sunbeam uses
 - `config.yml` contains configuration options that can be specified by the user when running an extension
 - `sbx_bracken.rules` contains the rules (logic/commands run) of the extension

 
    
