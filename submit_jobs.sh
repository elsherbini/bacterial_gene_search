#!/bin/bash

# Path to Snakefile
SNAKEFILE="workflow/Snakefile"

# Number of jobs to run in parallel
JOBS=10

# Profile to use (if any)
PROFILE=""

# Parse command line arguments
while getopts "j:p:" opt; do
  case $opt in
    j) JOBS=$OPTARG ;;
    p) PROFILE="--profile $OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Run snakemake
echo "Running snakemake with $JOBS jobs"
snakemake -s $SNAKEFILE --use-conda $PROFILE --jobs $JOBS 