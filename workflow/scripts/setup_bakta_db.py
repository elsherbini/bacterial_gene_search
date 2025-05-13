#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
from snakemake.shell import shell

# Get input variables from Snakemake
db_path = snakemake.params.db_path
download_db = snakemake.params.download_db

# Create directory if it doesn't exist
os.makedirs(db_path, exist_ok=True)

# Check if the database already exists (look for bakta.db at the specified path)
db_exists = os.path.exists(f"{db_path}/bakta.db")

if not db_exists:
    if download_db:
        print(f"Downloading Bakta database to {db_path}", file=sys.stderr)
        
        # Create a temporary directory for the download
        temp_download_dir = f"{db_path}_temp"
        os.makedirs(temp_download_dir, exist_ok=True)
        
        # Run the download command to the temporary directory
        shell(f"bakta_db download --output {temp_download_dir}")
        
        # Check if the download created the expected 'db' subdirectory
        db_subdir = f"{temp_download_dir}/db"
        if os.path.exists(db_subdir):
            # Move all files from the db subdirectory to the target directory
            print(f"Moving files from {db_subdir} to {db_path}", file=sys.stderr)
            for item in os.listdir(db_subdir):
                src = os.path.join(db_subdir, item)
                dst = os.path.join(db_path, item)
                shutil.move(src, dst)
            
            # Clean up the temporary directory
            shutil.rmtree(temp_download_dir)
        else:
            # If no db subdirectory was created, just move everything 
            print(f"No 'db' subdirectory found, moving all files from {temp_download_dir} to {db_path}", file=sys.stderr)
            for item in os.listdir(temp_download_dir):
                src = os.path.join(temp_download_dir, item)
                dst = os.path.join(db_path, item)
                shutil.move(src, dst)
            
            # Clean up the temporary directory
            shutil.rmtree(temp_download_dir)
        
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