#!/bin/bash

#SBATCH --account=pi-mstephens
#SBATCH --time=30:00:00
#SBATCH --ntasks=8                         # Run on 8 tasks
#SBATCH --cpus-per-task=8                  # Number of CPU cores (threads) per task
#SBATCH --mem=150G

# Get the input data file as specified by the command-line argument.
# INPUT=${1}

# For reproducibility, add version numbers to the "module load" calls.
module load R/4.4.1

# Run the R script.
echo "Running fash_analysis.r"
Rscript fash_analysis.R > fash_analysis.Rout