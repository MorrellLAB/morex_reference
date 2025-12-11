#!/usr/bin/env python3
"""
Title: MAF State Recovery (Complete Coverage)
Description: 
    Iterates through an AnchorWave MAF alignment to identify the nucleotide 
    state of a query genome (e.g., Hordeum bulbosum) at specific positions.
    
    Ensures ALL input VCF positions are reported.
    - If found in alignment: Reports base (A/C/G/T) or '-' (Deletion).
    - If NOT found in alignment: Reports 'N'.

Usage:
    python maf_state_recovery.py --maf alignment.maf --vcf variants.vcf.gz > hb_states.bed
"""

import sys
import gzip
import argparse
from Bio import AlignIO

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract query genome states from MAF to BED format."
    )
    parser.add_argument("--maf", required=True, help="Path to .maf alignment.")
    parser.add_argument("--vcf", required=True, help="Path to .vcf(.gz) file.")
    parser.add_argument("--output", default=sys.stdout, type=argparse.FileType('w'),
                        help="Output file path (default: stdout).")
    return parser.parse_args()

def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def load_vcf_targets(vcf_path):
    """
    Loads VCF positions into a dictionary.
    Returns: dict { 'chr1H': {100: {'ref': 'A', 'alt': 'G', 'id': 'rs1', 'found': False}} }
    """
    targets = {}
    with open_file(vcf_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError: continue
            
            snp_id = parts[2] if parts[2] else "."
            ref = parts[3]
            alt = parts[4]
            
            if chrom not in targets:
                targets[chrom] = {}
            
            # Use a dict to store details + a 'found' flag
            targets[chrom][pos] = {
                'ref': ref, 
                'alt': alt, 
                'id': snp_id, 
                'found': False
            }
    return targets

def process_alignment(maf_path, targets, output_handle):
    """
    Streams MAF, finds overlaps, prints them, and marks them as found.
    """
    for block in AlignIO.parse(maf_path, "maf"):
        ref_rec = block[0]
        qry_rec = block[1]
        
        # 1. Parse Chromosome
        full_ref_name = ref_rec.id
        ref_chrom = full_ref_name.split('.')[-1] if '.' in full_ref_name else full_ref_name

        if ref_chrom not in targets:
            continue

        # 2. Block Coordinates
        block_start = ref_rec.annotations['start'] 
        block_size = ref_rec.annotations['size']
        block_end = block_start + block_size

        # 3. Identify targets in block
        block_targets = []
        for t_pos in targets[ref_chrom]:
            # VCF is 1-based, MAF start is 0-based. 
            # Logic: block_start (0-based) < t_pos (1-based) <= block_end (0-based end)
            if (block_start + 1) <= t_pos <= block_end:
                # Check if already found (to handle duplicate blocks/overlaps if any)
                if not targets[ref_chrom][t_pos]['found']:
                    block_targets.append(t_pos)

        if not block_targets:
            continue
        
        block_targets.sort()
        
        # 4. Walk Alignment
        current_genomic_pos = block_start + 1 
        target_idx = 0
        seq_ref = str(ref_rec.seq)
        seq_qry = str(qry_rec.seq)
        
        for i, ref_base in enumerate(seq_ref):
            if target_idx >= len(block_targets):
                break
            
            target_pos = block_targets[target_idx]
            
            if current_genomic_pos == target_pos:
                hb_base = seq_qry[i]
                info = targets[ref_chrom][target_pos]
                
                # Mark as found so we don't print "N" later
                targets[ref_chrom][target_pos]['found'] = True
                
                # Determine state
                state = hb_base
                # You requested that 'DEL' not be used. 
                # Standard convention: '-' is a known deletion. 'N' is unknown.
                # If you want known deletions to be '-', keep this:
                if hb_base == '-':
                    state = '-' 
                
                out_line = [
                    ref_chrom, 
                    str(target_pos - 1), 
                    str(target_pos), 
                    info['id'],
                    info['ref'], 
                    info['alt'], 
                    state
                ]
                output_handle.write("\t".join(out_line) + "\n")
                target_idx += 1
            
            if ref_base != "-":
                current_genomic_pos += 1

def report_missing(targets, output_handle):
    """
    Iterates through all targets and prints 'N' for those not marked as found.
    """
    # To keep output somewhat sorted, we iterate explicitly
    sorted_chroms = sorted(targets.keys())
    
    for chrom in sorted_chroms:
        sorted_pos = sorted(targets[chrom].keys())
        for pos in sorted_pos:
            info = targets[chrom][pos]
            if not info['found']:
                # Not found in any MAF block -> N
                out_line = [
                    chrom,
                    str(pos - 1),
                    str(pos),
                    info['id'],
                    info['ref'],
                    info['alt'],
                    "N"
                ]
                output_handle.write("\t".join(out_line) + "\n")

def main():
    args = parse_arguments()
    
    sys.stderr.write(f"Loading targets from {args.vcf}...\n")
    targets = load_vcf_targets(args.vcf)
    
    sys.stderr.write(f"Scanning MAF {args.maf}...\n")
    process_alignment(args.maf, targets, args.output)
    
    sys.stderr.write("Reporting missing sites as 'N'...\n")
    report_missing(targets, args.output)
    
    sys.stderr.write("Done.\n")

if __name__ == "__main__":
    main()
