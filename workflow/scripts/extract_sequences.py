#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
import pandas as pd
import re

# Get variables from Snakemake
hmmsearch_output = snakemake.input.hmmsearch
protein_fasta = snakemake.input.proteins
gbff_file = snakemake.input.gbff
pfam_list_file = snakemake.input.pfam_list
output_fasta = snakemake.output.fasta
output_nucleotide = snakemake.output.nucleotide
gene_name = snakemake.wildcards.gene
genome_id = snakemake.wildcards.genome_id
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
        f.write(f"# No hits found for {gene_name} ({gene_description}) in {genome_id} with E-value <= {evalue_cutoff}\n")
    with open(output_nucleotide, 'w') as f:
        f.write(f"# No hits found for {gene_name} ({gene_description}) in {genome_id} with E-value <= {evalue_cutoff}\n")
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

# Build a mapping from protein IDs to nucleotide sequences
protein_to_nucleotide = {}
for record in SeqIO.parse(gbff_file, "genbank"):
    for feature in record.features:
        if feature.type == "CDS" and "protein_id" in feature.qualifiers:
            protein_id = feature.qualifiers["protein_id"][0]
            if protein_id in [hit_id for hit_id, _ in hits]:
                if feature.location:
                    start = feature.location.start
                    end = feature.location.end
                    strand = feature.location.strand
                    nuc_seq = feature.extract(record.seq)
                    protein_to_nucleotide[protein_id] = {
                        "seq": str(nuc_seq),
                        "location": f"{start}..{end}",
                        "strand": "+" if strand == 1 else "-",
                        "contig": record.id
                    }

# Extract nucleotide sequences
with open(output_nucleotide, 'w') as f:
    for hit_id, evalue in hits:
        if hit_id in protein_to_nucleotide:
            nuc_info = protein_to_nucleotide[hit_id]
            f.write(f">{hit_id} | {gene_name} ({gene_description}) | E-value={evalue:.2e} | {nuc_info['contig']}:{nuc_info['location']}({nuc_info['strand']})\n")
            f.write(f"{nuc_info['seq']}\n")
        else:
            f.write(f"# Could not find nucleotide sequence for protein {hit_id}\n")

print(f"Extracted {len(hits)} sequences for {gene_name} ({gene_description}) in {genome_id}", file=sys.stderr) 