# sbx_kraken2
This is an extension to the [Sunbeam pipeline](https://github.com/sunbeam-labs/sunbeam) to run the kraken2 read classifier. It works on qc'd reads, and requires that you have a kraken2 database available. Check out the [instructions here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#kraken-2-databases) for help with making a kraken2 database.

## Installing

Clone the repo into your sunbeam `extensions/` folder and add the new options to your existing configuration file. Make sure you've [installed Sunbeam](https://sunbeam.readthedocs.io/en/latest/quickstart.html) first!

    source activate sunbeam
    cd $SUNBEAM_DIR
    git clone https://github.com/louiejtaylor/sbx_kraken2/ extensions/sbx_kraken2

Add the options to your config file (replace "sunbeam_config.yml" with the name of your config file).

    cat extensions/sbx_kraken2/config.yml >> sunbeam_config.yml

## Configuration

In your configuration file, make sure you choose how many threads you'd like kraken2 to use, as well as providing the path to your database:

```
sbx_kraken2:
  threads: 4
  group_file: '/path/to/your/kraken2/database_dir'
```

## Running

Finally, run Sunbeam as usual (with the `--use-conda` flag) with your extension's target rule specified:

    sunbeam run --use-conda --configfile=sunbeam_config.yml all_kraken2

The `--use-conda` flag is crucial as it enables running kraken2 in an isolated environment.

This rule generates output in the following location: `sunbeam_output/classify/kraken2/`

## Contents

 - `sbx_kraken2_env.yml` specifies the extension's dependencies and provides the environment that Sunbeam uses
 - `config.yml` contains configuration options that can be specified by the user when running an extension
 - `sbx_kraken2.rules` contains the rules (logic/commands run) of the extension

 
    
 
