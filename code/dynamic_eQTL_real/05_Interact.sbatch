#!/bin/bash

#SBATCH --account=pi-mstephens
#SBATCH --time=36:00:00
#SBATCH --ntasks=1                         # Run on 1 task
#SBATCH --cpus-per-task=1                  # Number of CPU cores (threads) per task
#SBATCH --mem=150G

# Get the input data file as specified by the command-line argument.
# INPUT=${1}

# For reproducibility, add version numbers to the "module load" calls.
module load R/4.4.1

# Run the R script.
echo "Running 05_Interact.R"
Rscript 05_Interact.R > 05_Interact.Rout