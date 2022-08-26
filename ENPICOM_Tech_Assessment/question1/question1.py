#!/usr/bin/env python

# This script finds the 'best matching' sequence (b) out of a set of sequences (B), for a given sequence (a).
# The 'best matching' sequence is that which can be transformed into sequence (a) with the least amount of
# single character modifications (i.e. insertions, deletions, substitutions).

import argparse
import os.path
import numpy


def check_input_file(file):
    """ This function checks that the input file exists. """
    if not os.path.exists(file):
        return False
    return True


def find_best_match(a_file, b_set_dict):
    """
    Find the best matching sequence based on their Levenshtein distance to the target sequence (a).
    :param a_file: path to the file containing the target sequence (a).
    :param b_set_dict: dictionary of sequences (B) from where the best match (b) will be found.
    """
    # Open file containing 'a' in reading mode
    a = open(a_file, 'r')
    # Create a string to parse the lines corresponding to one sequence in the FASTA file
    a_sequence = ''
    for line in a:
        a_sequence += line.strip()
    # Create a dictionary to store the Levenshtein distance between 'a' and each sequence in 'B'
    levenshtein_distances = {}
    for sequence_id, b_sequence in b_set_dict.items():
        levenshtein_distances[sequence_id] = calculate_levenshtein_distance(a_sequence, b_sequence)
    # Find the best match by finding the sequence with the min Levenshtein distance to 'a'
    best_match = min(levenshtein_distances, key=levenshtein_distances.get)
    return best_match, str(levenshtein_distances[best_match]), b_set_dict[best_match]


def calculate_levenshtein_distance(a_sequence, b_sequence):
    """ This function calculates the Levenshtein distance between two sequences 'a' and 'b'."""
    # Create a null matrix for the calculations
    a_rows = len(a_sequence)
    b_cols = len(b_sequence)
    distances = numpy.zeros((a_rows + 1, b_cols + 1), dtype=int)
    # Fill in the first row and column of the matrix
    distances[:, 0] = [0] + list(range(1, a_rows + 1))
    distances[0, :] = [0] + list(range(1, b_cols + 1))
    # Calculate the Levenshtein distance
    for j in range(1, b_cols + 1):
        for i in range(1, a_rows + 1):
            if a_sequence[i-1] == b_sequence[j-1]:
                distance = 0
            else:
                distance = 1
            distances[i, j] = min(distances[i-1, j] + 1,
                                  distances[i, j-1] + 1,
                                  distances[i-1, j-1] + distance)
    return distances[i, j]


def create_b_set_dict(b_set_file):
    """ This function creates a dictionary of B with the IDs from the FASTA file as keys
     and their corresponding sequence as values."""
    # Open the database (B) FASTA file in reading mode
    b_set = open(b_set_file, 'r')
    b_set_dict = {}
    for line in b_set:
        line = line.strip()
        # If a line starts with '>' it's an ID and the start of a new sequence
        if line.startswith('>') and line not in b_set_dict:
            key = line.strip('>')
            b_set_dict[key] = ""
        # If the line does not start with '>' it's part of a sequence, so it's stored as a value
        else:
            b_set_dict[key] += line
    return b_set_dict


def parse_output_file(best_match):
    """ This function creates the output file."""
    with open('best_match.txt', 'w') as output:
        output.write(">" + best_match[0] + "\n")
        output.write("Single-character modifications: " + best_match[1].strip() + "\n")
        output.write(best_match[2] + "\n")


def main(a_file, b_set_file):
    # Check that the input files exist.
    if check_input_file(a_file) and check_input_file(b_set_file):
        # If the files exist, use the functions from above to find the best matching sequence (b).
        b_set_dict = create_b_set_dict(b_set_file)
        best_match = find_best_match(a_file, b_set_dict)
        parse_output_file(best_match)
    else:
        # If one of the files does not exist in the directory, ask the user to check.
        print(a_file + " or " + b_set_file + " does not exist in the directory. Please use existing FASTA files "
                                             "as input.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script finds the best matching sequence (b) to out of the set'
                                                 'of sequences (B), for a given sequence (a).'
                                                 'B and a are provided by the user '
                                                 'The output is a file containing a, b and the percentage of similarity'
                                                 'between a and b.')

    parser.add_argument("-a", "--target_sequence", help="Sequence for which the best matching sequence will be found",
                        required=True)
    parser.add_argument("-B", "--set_of_sequences", help="Set of sequences out of which the best matching pair will be"
                                                         "found.", required=True)
    args = parser.parse_args()

    a_filename = args.target_sequence
    B_filename = args.set_of_sequences

    """ Calls the "main" function.
    To make sure that it runs, the file must be run with the following command line:
        python3 question1.py -a a.fa -B B.fa"""
    main(a_filename, B_filename)
