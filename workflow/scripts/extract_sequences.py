#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
import pandas as pd
import re

# Get variables from Snakemake
hmmsearch_output = snakemake.input.hmmsearch
protein_fasta = snakemake.input.proteins
pfam_list_file = snakemake.input.pfam_list
output_fasta = snakemake.output.fasta
output_nucleotide = snakemake.output.get("nucleotide", None)
gene_name = snakemake.wildcards.gene
sample_id = snakemake.wildcards.sample
evalue_cutoff = float(snakemake.params.evalue_cutoff)

# Read Pfam descriptions
pfam_df = pd.read_csv(pfam_list_file)
gene_description = ""
for _, row in pfam_df.iterrows():
    if row['pfam_accession'] == gene_name:
        gene_description = row['description']
        break

# Parse hmmsearch output to get hits
hits = []
current_hit = None
with open(hmmsearch_output, 'r') as f:
    in_hits = False
    for line in f:
        if line.startswith(">>"):
            in_hits = True
            current_hit = line.strip().split()[1]
        elif in_hits and re.match(r'^\s+E-value', line):
            # Next line will have the E-value
            continue
        elif in_hits and current_hit and re.match(r'^\s+\d', line):
            # This is the score line with E-value
            fields = line.strip().split()
            evalue = float(fields[4])
            if evalue <= evalue_cutoff:
                hits.append((current_hit, evalue))
            current_hit = None

# If no hits found, create empty output files
if not hits:
    with open(output_fasta, 'w') as f:
        f.write(f"# No hits found for {gene_name} ({gene_description}) in {sample_id} with E-value <= {evalue_cutoff}\n")
    if output_nucleotide:
        with open(output_nucleotide, 'w') as f:
            f.write(f"# No hits found for {gene_name} ({gene_description}) in {sample_id} with E-value <= {evalue_cutoff}\n")
    sys.exit(0)

# Load protein sequences
proteins = SeqIO.to_dict(SeqIO.parse(protein_fasta, "fasta"))

# Extract hit sequences and write to output
with open(output_fasta, 'w') as f:
    for hit_id, evalue in hits:
        if hit_id in proteins:
            seq_record = proteins[hit_id]
            seq_record.description = f"{seq_record.description} | {gene_name} ({gene_description}) | E-value={evalue:.2e}"
            SeqIO.write(seq_record, f, "fasta")

# Extract nucleotide sequences if requested
if output_nucleotide:
    # We need to map from protein ID to nucleotide ID
    # This typically requires parsing the GBK file, but for this script we'll assume
    # the protein IDs can be mapped to nucleotide features in a GBK file
    # (In a real implementation, you would need to add this mapping logic)
    with open(output_nucleotide, 'w') as f:
        f.write(f"# Nucleotide sequences for {gene_name} ({gene_description}) in {sample_id}\n")
        f.write(f"# To properly implement this, we need to map protein IDs to nucleotide features\n")
        for hit_id, evalue in hits:
            f.write(f">{hit_id} | {gene_name} ({gene_description}) | E-value={evalue:.2e}\n")
            f.write("NUCLEOTIDE_SEQUENCE_PLACEHOLDER\n")

print(f"Extracted {len(hits)} sequences for {gene_name} ({gene_description}) in {sample_id}", file=sys.stderr) 