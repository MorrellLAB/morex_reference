#!/usr/bin/env python3

# After running SNP_Utils (https://github.com/mojaveazure/SNP_Utils), some
# SNPs will fail and not have any BLAST hits (sometimes because our identity
# threshold is very high). So, we use IPK BLAST server (with the latest morex
# reference genome) to identify the best possible hit and the associated
# physical positions.

# IMPORTANT CAVEAT: Currently, this snp uses contextual fasta sequences and
# assumes the SNP is in the middle (we have the same number of bases on
# both sides of the SNP)

import os
import sys
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from bs4 import BeautifulSoup

# IUPAC nucleotide base code: [Base, Reverse complement]
IUPAC_TABLE = {
    'R': [['A', 'G'], 'Y'],
    'Y': [['C', 'T'], 'R'],
    'S': [['G', 'C'], 'S'],
    'W': [['A', 'T'], 'W'],
    'K': [['T', 'G'], 'M'],
    'M': [['A', 'C'], 'K'],
    'B': [['C', 'G', 'T'], 'V'],
    'D': [['A', 'G', 'T'], 'H'],
    'H': [['A', 'C', 'T'], 'D'],
    'V': [['A', 'C', 'G'], 'B'],
    'N': [['A', 'C', 'G', 'T'], 'N']
}


def split(seq_str):
    """Split string."""
    temp = [char for char in seq_str]
    return temp


def read_fasta(fasta_file):
    """Parse the FASTA file and keep only the sequence identifiers and
    sequences in a dictionary."""
    fasta_dict = {}
    with open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_dict[record.id] = split(record.seq)
    return fasta_dict


def read_html(blast_results_html):
    """Read in html file containing IPK BLAST search results.
    Multiple SNP searches can be stored in a single HTML file."""
    with open(blast_results_html) as fp:
        soup = BeautifulSoup(''.join(fp), features="lxml")
    return soup


def closest(idx_lst, K):
    """Identify closest position to the middle of the sequence. For the BOPA and 9K
    SNPs, many of the SNPs are located in the middle of the contextual sequence, however
    this is not always the case. When the SNP is not in the center of the contextual
    sequence, pick the closest SNP."""
    return idx_lst[min(range(len(idx_lst)), key = lambda i: abs(idx_lst[i] - K))]


def context_snp_win(fasta_dict, snp_name):
    """Generate a 7bp window for the contextual sequence that contains the
    SNP. For BOPA and 9K SNPs, most of the time the SNP is in the middle
    of the contextual fasta sequence. If it is not in the middle, move
    left and right until we find the SNP."""
    # Identify the middle of the contextual sequence
    mid_context_snpidx = math.floor(len(fasta_dict[snp_name])/2)
    # Generate 7bp window containing SNP
    # This is case insensitive
    if fasta_dict[snp_name][mid_context_snpidx].upper() not in IUPAC_TABLE.keys():
        print("Middle of sequence is not the SNP, move left/right until we identify the SNP")
        # Move left from center until we find SNP
        for i in range(mid_context_snpidx, 0, -1):
            if fasta_dict[snp_name][i] not in IUPAC_TABLE.keys():
                continue
            else:
                # We found the SNP, save SNP index
                lsnp_idx = i
                break
        # Move right from the center until we find SNP
        for i in range(mid_context_snpidx, len(fasta_dict[snp_name])):
            if fasta_dict[snp_name][i].upper() not in IUPAC_TABLE.keys():
                continue
            else:
                # We found the SNP
                rsnp_idx = i
                break
        # Identify the position closest to center to use as SNP
        cand_snp_idx = [lsnp_idx, rsnp_idx]
        context_snpidx = closest(cand_snp_idx, mid_context_snpidx)
        context_seq = ''.join(fasta_dict[snp_name][context_snpidx-3: context_snpidx+4])
    else:
        # Center position is the SNP
        context_snpidx = mid_context_snpidx
        context_seq = ''.join(fasta_dict[snp_name][context_snpidx-3: context_snpidx+4])
    return (context_snpidx, context_seq)


def cat_query_seq(current_hit):
    """Takes in query sequence and concatenates sequence for
    easier search for SNP."""
    query_idx = []
    query_seq_list = []
    for i, elem in enumerate(current_hit):
        if elem.startswith('Query'):
            # Add index to query_idx
            query_idx.append(i)
            # Save query sequence
            query_seq_list.append(",".join(elem.split()).split(',')[2])
    # Concatenate the query sequence
    query_seq = "".join(query_seq_list)
    return query_seq


def pick_best_hit(current_snp, context_snp_seq):
    """Takes in a single SNP, picks the best BLAST hit that also contains
    the SNP in the contextual sequence, and returns data for the best
    BLAST hit for the current SNP."""
    # Step 1: create a list containing summaries of hits for current snp
    summary = []
    for counter, h in enumerate(current_snp):
        # This pattern depends on the html summary table and may need to
        # be modified accordingly
        if h.startswith('chr'):
            summary.append([counter, h.split()])
    # Pick the best hit
    # BLAST hits are always sorted by the best hit (highest bit score, lowest
    # E-value) first. We will use this assumption to pick the best SNP
    # IMPORTANT note to self: Double check this assumption!!!! And make sure
    # to link to documentation that specifies this!
    best_hit = summary[0]
    # Step 2: split hits for current snp into sublists for easier processing
    # Each SNP can have multiple "hits". "Hits" usually start with '>lcl'
    # and can have multiple matches (alignments)
    # Get indices for lines that start with '>lcl'
    hit_idx = []
    for counter, h in enumerate(current_snp):
        # This '>lcl' pattern will depend on the html file and may need to
        # be modified accordingly
        # Feature: Add check if no pattern is found, then exit with message
        #if h.startswith('>lcl'):
        if h.startswith('>chr'):
            hit_idx.append(counter)
    # Identify ending index of current SNP
    # Each SNP search ends with a few lines that starts with 'Lambda'
    tmp_end = []
    for i, elem in enumerate(current_snp):
        if "Lambda" in elem:
            tmp_end.append(i)
    # There are two 'Lambda' patterns that usually show up, we only need
    # the index of the first occurrence of the pattern
    snp_idx_end = tmp_end[0]
    # For current SNP, split into sublists based on lines starting with '>lcl'
    last_idx = len(hit_idx) - 1
    snp_hit_split = []
    for i in range(0, len(hit_idx)):
        if i < last_idx:
            current_hit = current_snp[hit_idx[i]:hit_idx[i+1]]
            snp_hit_split.append(current_hit)
        else:
            # If last index, that means this is the last match for this
            # SNP, use an ending index for the pattern 'Lambda'
            current_hit = current_snp[hit_idx[i]:snp_idx_end]
            snp_hit_split.append(current_hit)
    # Step 3: Split multiple matches for each "hit" and store in list of lists
    score_idx = []
    for i, elem in enumerate(snp_hit_split[0]):
        if elem.startswith(' Score'):
            score_idx.append(i)
    # For each chr, the best hits are also sorted by E-value by default
    # So, we can create a range of indices to process
    # We will pick the best hit that contains the SNP
    for i, elem in enumerate(score_idx):
        # If elem is not the last index
        if elem != score_idx[-1]:
            tmp_query_seq = cat_query_seq(snp_hit_split[0][score_idx[i]:score_idx[i+1]])
            print(tmp_query_seq)
            if context_snp_seq in tmp_query_seq:
                print("SNP is in the query sequence of this hit")
                score_start_idx = score_idx[i]
                score_end_idx = score_idx[i+1]
                break
        else:
            # We are at the last element or have only one hit
            # for this chromosome. Go until the end
            #print(snp_hit_split[0][score_idx[i]:])
            tmp_query_seq = cat_query_seq(snp_hit_split[0][score_idx[i]:])
            print(tmp_query_seq)
            if context_snp_seq in tmp_query_seq:
                print("SNP is in the query sequence of this hit")
                score_start_idx = score_idx[i]
                score_end_idx = len(snp_hit_split[0]) - 1
                break
    # Save chr and length info first
    best_hit = snp_hit_split[0][0:2]
    # Then add best hit containing SNP
    best_hit.extend(snp_hit_split[0][score_start_idx:score_end_idx])
    # if len(score_idx) > 1:
    #     best_hit.append(snp_hit_split[0][score_start_idx:score_end_idx])
    # else:
    #     # If we only have one hit for this chromosome,
    #     # # go until the end
    #     best_hit.append(snp_hit_split[0][score_idx[0]:])
    return best_hit


def plus_plus_strand(query_snp, ref_allele, iupac_table):
    if query_snp not in ['A', 'C', 'T', 'G']:
        if query_snp in ['B', 'D', 'H', 'V']:
            print("We have more than 1 alternate alleles")
            # For alt allele, use the one that is not the ref allele
            alt_allele = []
            for i, elem in enumerate(iupac_table[query_snp][0]):
                # We only consider 2 nucleotides case
                if elem != ref_allele:
                    alt_allele.append(elem)
        else:
            # For alt allele, use the one that is not the ref allele
            for i, elem in enumerate(iupac_table[query_snp][0]):
                # We only consider 2 nucleotides case
                if elem != ref_allele:
                    alt_allele = elem
    return alt_allele


def plus_minus_strand(query_snp, ref_allele, iupac_table):
    if query_snp not in ['A', 'C', 'T', 'G']:
        if query_snp in ['B', 'D', 'H', 'V']:
            print("We have more than 1 alternate alleles")
            # For alt allele, use the one that is not the ref allele
            alt_allele = []
            for i, elem in enumerate(iupac_table[query_snp][0]):
                # We only consider 2 nucleotides case
                if elem != ref_allele:
                    alt_allele.append(elem)
            print("Ref allele:", ref_allele)
            print("Alt allele:", alt_allele)
        else:
            # For alt allele, use the one that is not the ref allele
            for i, elem in enumerate(iupac_table[query_snp][0]):
                # We only consider 2 nucleotides case
                if elem != ref_allele:
                    alt_allele = elem
            print("Ref allele:", ref_allele)
            print("Alt allele:", alt_allele)
    # Now, take the reverse complement
    rc_ref_allele = Seq(ref_allele, generic_dna).reverse_complement()[0]
    print("Rev comp ref allele:", rc_ref_allele)
    rc_alt_allele = Seq(alt_allele, generic_dna).reverse_complement()[0]
    print("Rev comp alt allele:", rc_alt_allele)
    return (rc_ref_allele, rc_alt_allele)


def extract_info(snp_name, current_best_hit, fasta_dict):
    """ """
    # Save chromosome info
    chrom = current_best_hit[0].strip('>').strip()
    # Save percent identity info
    #identity = current_best_hit[3].strip().split(',')[0].split(' ')[3]
    # Let's get a list of indices in best hit where the line start with
    # "Query" or "Sbjct"
    query_idx = []
    sbjct_idx = []
    query_seq_list = []
    sbjct_seq_list = []
    for i, elem in enumerate(current_best_hit):
        # Let's store info about the % identity and strand
        if elem.startswith(' Identities'):
            tmp_idt = elem
        if elem.startswith(' Strand'):
            tmp_strand = elem
        if elem.startswith('Query'):
            # Add index to query_idx
            query_idx.append(i)
            # Save query sequence
            query_seq_list.append(",".join(elem.split()).split(',')[2])
        if elem.startswith('Sbjct'):
            # Add index to sbjct_idx
            sbjct_idx.append(i)
            # Save subject sequence
            sbjct_seq_list.append(",".join(elem.split()).split(',')[2])
    # Pull out relevant parts and save for later
    chrom = current_best_hit[0].strip('>').strip()
    idt = tmp_idt.strip().split(',')[0].replace(" ", "")
    strand = tmp_strand.split('=')[1]
    info_field = "".join(["B;", idt, ",failed"])
    qual = str('.')
    filter_field = str('.')
    # Concatenate the sequences
    query_seq = "".join(query_seq_list)
    sbjct_seq = "".join(sbjct_seq_list)
    # Keep starting position of subject (ref)
    sbjct_start_pos = ",".join(current_best_hit[sbjct_idx[0]].split()).split(',')[1]
    # Find the physical position of the SNP in the reference
    # This assumes Plus/Plus strand
    context_seq = fasta_dict[snp_name][1]
    if strand == "Plus/Plus":
        if query_seq.find(context_seq) != -1:
            # Return the leftmost index of the context_seq + floor(len(context_seq)/2)
            qsnp_idx = query_seq.find(context_seq) + math.floor(len(context_seq)/2)
            query_snp = query_seq[qsnp_idx]
            # Now, we have the index for the SNP
            # Let's get the associated position in the reference (Sbjct)
            ref_allele = sbjct_seq[qsnp_idx]
            if ref_allele == "-":
                print("Reference has an insertion, manually fix the position for SNP:", snp_name, "\n")
            # Count number of indels that occur prior to SNP
            num_indels = sbjct_seq[:qsnp_idx].count('-')
            # To get the correct reference position, we need to subtract the number of indels
            # that occur prior to the SNP
            subject_pos = int(sbjct_start_pos) + qsnp_idx - num_indels
            # Identify alternate allele
            alt_allele = plus_plus_strand(query_snp, ref_allele, IUPAC_TABLE)
        else:
            print("Could not resolve position for SNP", snp_name, ", saving to log file\n")
            # Add code here to save SNP to log file
    elif strand == "Plus/Minus":
        print("coordinates are reversed")
        if query_seq.find(context_seq) != -1:
            # Return the leftmost index of the context_seq + floor(len(context_seq)/2)
            qsnp_idx = query_seq.find(context_seq) + math.floor(len(context_seq)/2)
            query_snp = query_seq[qsnp_idx]
            # Now, we have the index for the SNP
            # Let's get the associated position in the reference (Sbjct)
            ref_allele = sbjct_seq[qsnp_idx]
            print("Ref", ref_allele)
            # Count number of indels that occur prior to SNP
            num_indels = sbjct_seq[:qsnp_idx].count('-')
            # To get the correct reference position, we need to add the number of indels
            # due to Plus/Minus strand
            subject_pos = int(sbjct_start_pos) - qsnp_idx + num_indels
            # Identify alternate allele and take reverse complement
            ref_allele, alt_allele = plus_minus_strand(query_snp, ref_allele, IUPAC_TABLE)
            print("RC Ref", ref_allele)
            print("RC Alt", alt_allele)
    else:
        print("Strand is Minus/Plus, saving to log file...")
    # Save VCF line
    return [chrom, str(subject_pos), snp_name, ref_allele, alt_allele, qual, filter_field, info_field]


def main(FASTA, BLAST_RESULTS_HTML, OUT_FILE):
    """Driver function."""
    # Read in fasta file
    fasta_dict = read_fasta(os.path.expanduser(FASTA))
    # Read in HTML file that contains IPK BLAST search results
    # where the HTML contains search results for one or more SNP searches
    soup = read_html(os.path.expanduser(BLAST_RESULTS_HTML))
    # Process soup object and get it into a workable data structure.
    # Extract text
    text = soup.get_text()
    # Split by newline delimiter
    content = text.split('\n')
    # Identify start and end indices of SNPs (elements starting with 'Query=')
    start_idx = [i for i, j in enumerate(content) if j.startswith('Query=')]
    end_idx = [i for i, j in enumerate(content) if j.startswith('Effective search space used:')]
    # Create list of lists for each SNP
    l = []
    for i in range(0, len(start_idx)):
        tmp = content[start_idx[i]:end_idx[i]+1]
        l.append(tmp)
    # Generate 7bp windows of contextual fasta containing the SNP
    fasta_win_dict = {}
    for i, elem in enumerate(fasta_dict.keys()):
        tmp_context_snpidx, tmp_context_seq = context_snp_win(fasta_dict, elem)
        fasta_win_dict[elem] = [tmp_context_snpidx, tmp_context_seq]
    # Pick the best hit for each SNP using the pick_best_hit function
    bhs = {}
    for i, elem in enumerate(l):
        # Save current snp name
        csnp_name = elem[0].split()[1]
        # Probably need to store the output from pick_best_hist somehow
        cbhs = pick_best_hit(current_snp=elem, context_snp_seq=fasta_win_dict[csnp_name][1])
        # Add new dictionary key,value pair
        bhs[csnp_name] = cbhs
    # Start from clean file, check if file exists
    if os.path.exists(os.path.expanduser(OUT_FILE)):
        os.remove(os.path.expanduser(OUT_FILE))
    # Extract info from best hit and save to file
    with open(os.path.expanduser(OUT_FILE), 'a') as f:
        # Add header line
        f.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]) + "\n")
        for i, elem in enumerate(bhs):
            vcf_line = extract_info(snp_name=elem, current_best_hit=bhs[elem], fasta_dict=fasta_win_dict)
            # Save VCF line to file
            f.write("\t".join(vcf_line) + "\n")
    return

main(sys.argv[1], sys.argv[2], sys.argv[3])  # Run the program
