#!/usr/bin/env python3
"""
Script to rename FASTQ files to a standardized format:
- test_R1_001.fastq.gz -> test_001_R1.fastq.gz
- test_R2_001.fastq.gz -> test_001_R2.fastq.gz
- test.fastq.gz -> test_R1.fastq.gz
"""

import os
import re
import glob
import argparse

def rename_fastq_files(data_dir, dry_run=False):
    """
    Rename FASTQ files in the specified directory.
    
    Args:
        data_dir: Directory containing FASTQ files
        dry_run: If True, only print what would be renamed without actually renaming
    """
    files = glob.glob(os.path.join(data_dir, "*.fastq.gz"))
    renamed_count = 0
    
    for filepath in files:
        filename = os.path.basename(filepath)
        new_filename = None
        
        # Pattern 1: name_R1_suffix.fastq.gz or name_R2_suffix.fastq.gz
        match = re.search(r'^(.+)_R([12])_(.+)(\.fastq\.gz)$', filename)
        if match:
            base_name = match.group(1)
            read_num = match.group(2)
            suffix = match.group(3)
            extension = match.group(4)
            new_filename = f"{base_name}_{suffix}_R{read_num}{extension}"
        
        # Pattern 2: name.fastq.gz (no R1/R2 indicator, assume R1)
        elif re.search(r'^(.+)(\.fastq\.gz)$', filename) and '_R1' not in filename and '_R2' not in filename:
            match = re.search(r'^(.+)(\.fastq\.gz)$', filename)
            base_name = match.group(1)
            extension = match.group(2)
            new_filename = f"{base_name}_R1{extension}"
        
        # Perform rename if needed
        if new_filename and new_filename != filename:
            old_path = filepath
            new_path = os.path.join(data_dir, new_filename)
            
            if dry_run:
                print(f"[DRY RUN] Would rename: {filename} -> {new_filename}")
                renamed_count += 1
            else:
                if os.path.exists(new_path):
                    print(f"[SKIP] Target already exists: {new_filename}")
                else:
                    os.rename(old_path, new_path)
                    print(f"[RENAMED] {filename} -> {new_filename}")
                    renamed_count += 1
    
    if dry_run:
        print(f"\n[DRY RUN] Would rename {renamed_count} files")
    else:
        print(f"\n[DONE] Renamed {renamed_count} files")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename FASTQ files to standardized format")
    parser.add_argument("-d", "--directory", default="fastq/", 
                        help="Directory containing FASTQ files (default: fastq/)")
    parser.add_argument("-n", "--dry-run", action="store_true",
                        help="Show what would be renamed without actually renaming")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.directory):
        print(f"Error: Directory '{args.directory}' does not exist")
        exit(1)
    
    rename_fastq_files(args.directory, args.dry_run)