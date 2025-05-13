#!/usr/bin/env python3

import os
import sys
import subprocess
import tempfile
import re
import pandas as pd
from snakemake.shell import shell

def find_full_accession(base_acc, hmm_db):
    """Find the full accession with version number in the HMM database"""
    try:
        # Use grep to find the accession line
        cmd = f"grep -m 1 '^ACC {base_acc}\\.' {hmm_db}"
        result = subprocess.run(cmd, shell=True, check=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                               text=True)
        
        # Extract the full accession with version number
        print(result.stdout)
        match = re.search(r'ACC\s+(\S+)', result.stdout)
        if match:
            return match.group(1)
        return None
    except subprocess.CalledProcessError:
        return None

def smart_hmmfetch(hmm_db, accession, output_file):
    """Run hmmfetch with version-aware accession handling"""
    try:
        # First try with the exact accession
        subprocess.run(["hmmfetch", hmm_db, accession], 
                      stdout=open(output_file, 'w'), 
                      stderr=subprocess.PIPE, 
                      check=True)
        return True
    except subprocess.CalledProcessError:
        # If exact match fails, try to find the versioned accession
        full_acc = find_full_accession(accession, hmm_db)
        
        if full_acc:
            print(f"Found versioned accession: {full_acc} for {accession}", file=sys.stderr)
            try:
                subprocess.run(["hmmfetch", hmm_db, full_acc], 
                              stdout=open(output_file, 'w'), 
                              stderr=subprocess.PIPE, 
                              check=True)
                return True
            except subprocess.CalledProcessError:
                print(f"Error running hmmfetch with accession {full_acc}", file=sys.stderr)
                return False
        else:
            print(f"Accession {accession} not found in database", file=sys.stderr)
            return False

# Get variables from Snakemake
pfam_db = snakemake.params.pfam_db
pfam_list_file = snakemake.input.pfam_list
output_hmm = snakemake.output.hmm

# Read Pfam accession list from CSV
pfam_df = pd.read_csv(pfam_list_file)
pfam_accessions = pfam_df['pfam_accession'].tolist()

# Create a temporary directory to store individual HMM files
with tempfile.TemporaryDirectory() as temp_dir:
    print(f"Extracting {len(pfam_accessions)} HMMs from Pfam database", file=sys.stderr)
    
    # Extract each HMM individually using smart_hmmfetch
    successful = 0
    for i, acc in enumerate(pfam_accessions):
        temp_hmm = os.path.join(temp_dir, f"{acc}.hmm")
        if smart_hmmfetch(f"{pfam_db}/Pfam-A.hmm", acc, temp_hmm):
            successful += 1
    
    # Combine all successfully extracted HMMs into one file
    with open(output_hmm, 'w') as out_file:
        for acc in pfam_accessions:
            temp_hmm = os.path.join(temp_dir, f"{acc}.hmm")
            if os.path.exists(temp_hmm) and os.path.getsize(temp_hmm) > 0:
                with open(temp_hmm, 'r') as in_file:
                    out_file.write(in_file.read())
                    out_file.write("\n")  # Add newline between HMMs

print(f"Successfully extracted {successful}/{len(pfam_accessions)} HMMs to {output_hmm}", file=sys.stderr) 