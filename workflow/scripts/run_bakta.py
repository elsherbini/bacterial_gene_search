#!/usr/bin/env python3

import os
import shutil
import subprocess
import glob
from snakemake.shell import shell

# Get variables from Snakemake
genome_path = snakemake.input.genome
genome_type = snakemake.params.genome_type
db_path = snakemake.params.db_path
prefix = snakemake.params.prefix
outdir = snakemake.params.outdir
threads = snakemake.threads
output_gbff = snakemake.output.gbff
output_faa = snakemake.output.faa

# Create output directory if it doesn't exist
os.makedirs(outdir, exist_ok=True)

# Handle different genome types
if genome_type == 'bakta_dir':
    # For already Bakta-annotated genomes, copy or link the necessary files
    input_dir = genome_path
    
    # Find and copy the .gbff file
    gbff_files = glob.glob(os.path.join(input_dir, "*.gbff"))
    if gbff_files:
        shutil.copy(gbff_files[0], output_gbff)
    else:
        raise ValueError(f"No .gbff file found in Bakta directory: {input_dir}")
    
    # Find and copy the .faa file
    faa_files = glob.glob(os.path.join(input_dir, "*.faa"))
    if faa_files:
        shutil.copy(faa_files[0], output_faa)
    else:
        raise ValueError(f"No .faa file found in Bakta directory: {input_dir}")
else:
    # For unannotated genomes (fna), run bakta for annotation
    cmd = [
        "bakta",
        "--db", db_path,
        "--output", outdir,
        "--prefix", prefix,
        "--keep-contig-headers",
        genome_path,
        "--threads", str(threads)
    ]
    
    subprocess.run(cmd, check=True) 