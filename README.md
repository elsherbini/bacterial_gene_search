# Bacterial Gene Search

This workflow is designed for taking a list of one or more genomes and a list of one or more genes to search for. It returns a fasta file for each gene, and an alignment, nucleotide tree, and amino acid tree for each gene if you had more than one genome searched.

This tool uses bakta for genome annotation, the hmmer suite for gene search, mafft for gene alignment, and fasttree for making gene trees.

## Overview

The workflow performs the following steps:
1. Sets up a local Pfam database (if not already available)
2. Extracts the HMMs for the target genes from the Pfam database
3. Annotates input genomes with Bakta (if not already annotated)
4. Searches for the target genes in the annotated genomes using HMMER
5. Extracts the identified gene sequences (both amino acid and nucleotide)
6. For multiple genomes, creates alignments and phylogenetic trees

## Input Files

### Genome List

The input genomes are specified in a CSV file with the following columns:
- `genome_id`: Unique genome identifier
- `genome_path`: Path to the genome file (.fna assembly file) or directory containing Bakta results
- `genome_type`: Type of the genome file (`fna` or `bakta_dir`)

Example:
```
genome_id,genome_path,genome_type
sample1,resources/genomes/sample1.fna,fna
sample2,resources/bakta_results/sample2,bakta_dir
```

### Gene List

The input genes are specified as a CSV file with Pfam accessions and descriptions:
```
pfam_accession,description
PF00005,ABC transporter
PF00072,Response regulator receiver domain
PF00106,Short chain dehydrogenase
```

The `pfam_accession` column is required and must contain valid Pfam IDs (e.g., PF00005). The `description` column is optional but recommended for better labeling of output files and sequences.

## Configuration

The configuration is stored in `config/config.yaml` and includes the following sections:

- Input/output paths
- Database configurations
- Tool parameters (e.g., HMMER E-value cutoff)
- Resource configurations (threads, memory, etc.)

## Running the Pipeline

### Prerequisites

1. Snakemake (version ≥ 8.20.0)
2. Conda or Mamba package manager

### Execution

1. Clone this repository:
```
git clone <repository-url>
cd bacterial_gene_search
```

2. Edit the config files:
   - `config/config.yaml`: Set paths and parameters
   - `config/genome_list.csv`: Add your genomes
   - `config/pfam_list.csv`: Add Pfam accessions and descriptions for target genes

3. Run the pipeline:
```
./submit_jobs.sh
```

Or with specific parameters:
```
./submit_jobs.sh -j 20 -p slurm
```

### Options for submit_jobs.sh

- `-j <int>`: Number of jobs to run in parallel (default: 10)
- `-p <str>`: Snakemake profile to use (optional)

## Output Files

The pipeline generates the following output files:

- `results/sequences/{sample}/{gene}.faa`: Amino acid sequences for each gene in each sample
- `results/sequences/{sample}/{gene}.fna`: Nucleotide sequences for each gene in each sample
- `results/combined/{gene}.faa`: Combined amino acid sequences for each gene across all samples
- `results/combined/{gene}.fna`: Combined nucleotide sequences for each gene across all samples
- `results/alignments/{gene}.aln`: Multiple sequence alignment for each gene
- `results/trees/{gene}.aa.tree`: Amino acid-based phylogenetic tree for each gene
- `results/trees/{gene}.nuc.tree`: Nucleotide-based phylogenetic tree for each gene

## Notes

- The E-value cutoff for HMMER searches is configurable (default: 1e-10)
- For Bakta annotation, a local database is required. See the [Bakta documentation](https://github.com/oschwengers/bakta) for details on setting up the database.
- The workflow will automatically download and set up the Pfam database if it's not already available at the specified location.
