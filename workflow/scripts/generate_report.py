#!/usr/bin/env python3

import os
import sys
import pandas as pd
import json
from pathlib import Path

# Get variables from Snakemake
pfam_list = snakemake.input.pfam_list
results_dir = snakemake.params.results_dir
output_json = snakemake.output.json

# Read Pfam info
pfam_df = pd.read_csv(pfam_list)
pfams = pfam_df.to_dict(orient='records')

# Get samples
samples = snakemake.params.samples

# Initialize results dictionary
report_data = {
    "pfams": pfams,
    "samples": samples,
    "sample_count": len(samples),
    "genes": []
}

# For each gene, check presence in each sample
for pfam in pfams:
    pfam_id = pfam['pfam_accession']
    
    # Check which samples have this gene
    found_in = []
    for sample in samples:
        faa_file = os.path.join(results_dir, "sequences", sample, f"{pfam_id}.faa")
        
        # Check if file exists and is not empty (or just has a comment)
        if os.path.exists(faa_file):
            with open(faa_file, 'r') as f:
                content = f.read().strip()
                if content and not all(line.startswith('#') for line in content.split('\n')):
                    found_in.append(sample)
    
    # Check if tree is available
    tree_file = os.path.join(results_dir, "trees", f"{pfam_id}.aa.tree")
    tree_available = False
    if os.path.exists(tree_file):
        with open(tree_file, 'r') as f:
            content = f.read().strip()
            if content and not content.startswith('#'):
                tree_available = True
    
    # Add gene info to report data
    gene_info = {
        "pfam_accession": pfam_id,
        "description": pfam['description'],
        "found_in": found_in,
        "tree_available": tree_available
    }
    report_data["genes"].append(gene_info)

# Write report data to JSON file
with open(output_json, 'w') as f:
    json.dump(report_data, f, indent=2) 