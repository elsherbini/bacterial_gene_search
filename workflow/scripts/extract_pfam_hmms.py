#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
from snakemake.shell import shell

# Get variables from Snakemake
pfam_db = snakemake.params.pfam_db
pfam_list_file = snakemake.input.pfam_list
output_hmm = snakemake.output.hmm

# Read Pfam accession list from CSV
pfam_df = pd.read_csv(pfam_list_file)
pfam_accessions = pfam_df['pfam_accession'].tolist()

# Create temporary file with accessions
temp_file = snakemake.output.hmm + ".accessions.tmp"
with open(temp_file, 'w') as f:
    for acc in pfam_accessions:
        f.write(f"{acc}\n")

# Use hmmfetch to extract HMMs
shell(f"hmmfetch -f {pfam_db}/Pfam-A.hmm {temp_file} > {output_hmm}")

# Clean up temp file
os.remove(temp_file)

print(f"Extracted {len(pfam_accessions)} HMMs to {output_hmm}", file=sys.stderr) 