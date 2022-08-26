#!/usr/bin/env/ python3

import sys


def parse_fasta(fasta_file_contents: str):
    """ Parses a given fasta file, to extract sequences and their headers

    Keyword arguments:
    fasta_file_contents -- String of the content of a FASTA file

    Returns:
    Dictionary with headers as keys, and sequences as values
    """
    fasta_dict = {}
    header = ''  # Init header in function scope
    lines = fasta_file_contents.split('\n')
    for line in lines:
        if line.startswith('>'):
            header = line[1:]
            if header not in fasta_dict:
                fasta_dict[header] = ''
        else:
            # It's a sequence
            fasta_dict[header] += line.strip()
    return fasta_dict


def rev_comp(dna: str):
    """ Generates the reverse complement of a uppercase DNA string.

    Keyword arguments:
    dna -- String of the DNA

    Returns:
    String of the respective reverse complement
    """
    dna_base_pairs = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    reverse_complement = [dna_base_pairs[base] for base in reversed(dna)]
    return ''.join(reverse_complement)


def split_correct_reads(sequences: list):
    rev_comp_sequences = [rev_comp(x) for x in sequences]

    correct_reads = [x for x in sequences if sequences.count(x) + rev_comp_sequences.count(x) > 1]
    incorrect_reads = [x for x in sequences if x not in correct_reads]

    return correct_reads, incorrect_reads


def calc_hamm_dist(read_one: str, read_two: str):
    """ Calculates the hamming distance, i.e. the number of different bases,
    for 2 DNA sequences.
    Reads are assumed to be of equal length

    Keyword arguments:
    read_one -- First DNA sequence
    read_two -- Second DNA sequence, that is compared to the first one

    Returns:
    int, number of different bases
    """
    num_diff = 0
    for base_one, base_two in zip(read_one, read_two):
        if base_one != base_two:
            num_diff += 1
    return num_diff


def find_corrections(incorrect_reads: list, correct_reads: list):
    corrected_reads = []
    for read_forw in incorrect_reads:
        # Find the correct read
        for correct_read in correct_reads:
            if calc_hamm_dist(read_forw, correct_read) == 1:
                # Found it, add a tuple to the stored list
                corrected_reads.append((read_forw, correct_read))
                break
            elif calc_hamm_dist(rev_comp(read_forw), correct_read) == 1:
                # Found it, add a tuple to the stored list
                corrected_reads.append((read_forw, rev_comp(correct_read)))
                break
    return corrected_reads


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit('Usage: python3 read_correction.py /path/to/file.fasta')
    # First parse the FASTA file, to extract the sequences
    fasta_content = open(sys.argv[1]).read()
    fasta_dict = parse_fasta(fasta_content)
    forward_sequences = list(fasta_dict.values())  # Discard the headers for this exercise

    # Split out the duplicates, and thus the correct reads
    correct_reads, incorrect_reads = split_correct_reads(forward_sequences)
    # Now correct the incorrect reads
    corrections = find_corrections(incorrect_reads, correct_reads)
    # Last, print the corrections
    for original, correction in corrections:
        print("{}->{}".format(original, correction))
