#!/bin/bash
#SBATCH --job-name=DatabaseBuilding
#SBATCH --account=of33
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=50000
#SBATCH --time=5-00:00:00

# Load blast
module load blast

# Download taxonomy
../KRAKEN2/kraken2-build --download-taxonomy --db Kraken20201214

# Download libraries
../KRAKEN2/kraken2-build --download-library archaea --db Kraken20201214
../KRAKEN2/kraken2-build --download-library bacteria --db Kraken20201214
../KRAKEN2/kraken2-build --download-library fungi --db Kraken20201214
../KRAKEN2/kraken2-build --download-library protozoa --db Kraken20201214
../KRAKEN2/kraken2-build --download-library viral --db Kraken20201214

# Build database
../KRAKEN2/kraken2-build --build --db Kraken20201214
