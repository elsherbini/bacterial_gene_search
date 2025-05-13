#!/usr/bin/env python3

import argparse
import subprocess
import re
import sys

def find_full_accession(base_acc, hmm_db):
    """Find the full accession with version number in the HMM database"""
    try:
        # Use grep to find the accession line
        cmd = f"grep -m 1 '^ACC {base_acc}\\.' {hmm_db}"
        result = subprocess.run(cmd, shell=True, check=True, 
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                               text=True)
        
        # Extract the full accession with version number
        match = re.search(r'ACC\s+(\S+)', result.stdout)
        if match:
            return match.group(1)
        return None
    except subprocess.CalledProcessError:
        return None

def run_hmmfetch(hmm_db, accession, output=None):
    """Run hmmfetch with the given accession"""
    cmd = ["hmmfetch", hmm_db, accession]
    
    if output:
        with open(output, 'w') as out_file:
            subprocess.run(cmd, stdout=out_file, check=True)
    else:
        subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description='Smart hmmfetch that handles Pfam accession versions')
    parser.add_argument('hmm_db', help='Path to the HMM database')
    parser.add_argument('accession', help='Base accession without version (e.g. PF09867)')
    parser.add_argument('-o', '--output', help='Output file (optional)')
    
    args = parser.parse_args()
    
    # First try with the exact accession
    try:
        if args.output:
            run_hmmfetch(args.hmm_db, args.accession, args.output)
        else:
            run_hmmfetch(args.hmm_db, args.accession)
        return 0
    except subprocess.CalledProcessError:
        # If exact match fails, try to find the versioned accession
        full_acc = find_full_accession(args.accession, args.hmm_db)
        
        if full_acc:
            print(f"Found versioned accession: {full_acc}", file=sys.stderr)
            try:
                if args.output:
                    run_hmmfetch(args.hmm_db, full_acc, args.output)
                else:
                    run_hmmfetch(args.hmm_db, full_acc)
                return 0
            except subprocess.CalledProcessError:
                print(f"Error running hmmfetch with accession {full_acc}", file=sys.stderr)
                return 1
        else:
            print(f"Accession {args.accession} not found in database", file=sys.stderr)
            return 1

if __name__ == "__main__":
    sys.exit(main()) 