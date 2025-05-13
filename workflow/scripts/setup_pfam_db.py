#!/usr/bin/env python3

import os
import subprocess
import sys
from snakemake.shell import shell

# Get input variables from Snakemake
pfam_db_path = snakemake.output[0].rstrip(".snakemake_setup_complete")
pfam_download_url = snakemake.params.download_url

# Create directory if it doesn't exist
os.makedirs(os.path.dirname(pfam_db_path), exist_ok=True)

# Download Pfam database if it doesn't exist
if not os.path.exists(pfam_db_path):
    print(f"Downloading Pfam database from {pfam_download_url}", file=sys.stderr)
    
    # Download gzipped file
    shell(f"wget -O {pfam_db_path}.gz {pfam_download_url}")
    
    # Decompress
    shell(f"gzip -d {pfam_db_path}.gz")
    
    # Press database for faster searches
    shell(f"hmmpress {pfam_db_path}")
else:
    print(f"Pfam database already exists at {pfam_db_path}", file=sys.stderr)
    
    # Check if pressed files exist, if not press again
    if not all(os.path.exists(f"{pfam_db_path}.{ext}") for ext in ["h3f", "h3i", "h3m", "h3p"]):
        print("Re-pressing Pfam database", file=sys.stderr)
        shell(f"hmmpress {pfam_db_path}")

print("Pfam database setup complete", file=sys.stderr) 