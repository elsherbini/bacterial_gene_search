#!/usr/bin/env python3

import os
import sys
import subprocess
from snakemake.shell import shell

# Get input variables from Snakemake
db_path = snakemake.params.db_path
download_db = snakemake.params.download_db

# Create directory if it doesn't exist
os.makedirs(os.path.dirname(db_path), exist_ok=True)

# Check if the database already exists
db_exists = os.path.exists(f"{db_path}/db.json")

if not db_exists:
    if download_db:
        print(f"Downloading Bakta database to {db_path}", file=sys.stderr)
        shell(f"bakta_db download --output {db_path}")
        
        # Verify download succeeded
        if not os.path.exists(f"{db_path}/db.json"):
            sys.exit(f"ERROR: Failed to download Bakta database to {db_path}")
    else:
        sys.exit(f"ERROR: Bakta database not found at {db_path} and download_db is set to false. "
                 f"Please either download the database manually or set download_db to true in the config.")
else:
    print(f"Bakta database already exists at {db_path}", file=sys.stderr)

# Create a flag file to mark successful setup
with open(snakemake.output[0], "w") as f:
    f.write("Database setup completed at " + db_path + "\n")

print("Bakta database setup complete", file=sys.stderr) 