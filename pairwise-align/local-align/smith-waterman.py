""" 
LocallSequenceAlignment

Summary
-------
Performs local DNA sequence alignment using a similarity scoring scheme.
""" 

import sys
import numpy as np

def parse_FASTA_file(input):
    """ Parses FASTA files containing 2 nucleotide sequences 
    """

    s1 = ""
    s2 = ""
    s = ""

    for line in input:
        if line.startswith(">"):  
            s += "!" 
        else:
            stripped = line.rstrip("\n")
            s += stripped

    s_list = s.split("!") 
    s1 = s_list[1] 
    s2 = s_list[2]

    return(s1, s2)


def local_align(seq1, seq2, M, m, g):
    """ Implementation of the Smith-Waterman Algorithm 

    Example Usage
    -------------
    # Given a defined scoring scheme, get optimal global sequence alignments
    alignments = local_align(seq1, seq2, match, mismatch, gap)

    Params
    -------
    seq1: sequence 1
    seq2: sequence 2
    M: match score
    m: mismatch penalty
    g: gap penalty
    """

    # Create a matrix to store the alignment scores
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = [[0] * cols for _ in range(rows)]

    # Variables to store maximum score and its position in the matrix
    max_score = 0
    max_pos = []

    # Fill the score matrix
    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i - 1] == seq2[j - 1]:
                match = score_matrix[i - 1][j - 1] + M
            else:
                match = score_matrix[i - 1][j - 1] + m
            delete = score_matrix[i - 1][j] + g
            insert = score_matrix[i][j - 1] + g

            # Take the maximum of match, delete, insert, and 0 (local alignment)
            score_matrix[i][j] = max(match, delete, insert, 0)

            # Update the maximum score and its poss
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = [(i, j)]
            elif score_matrix[i][j] == max_score:
                max_pos.append((i, j))

    # Trace back to find the local alignment strings
    alignments = []
    for pos in max_pos:
        optimal_seq1 = ""
        optimal_seq2 = ""
        i, j = pos
        while score_matrix[i][j] != 0:
            if seq1[i - 1] == seq2[j - 1]:
                optimal_seq1 = seq1[i - 1] + optimal_seq1
                optimal_seq2 = seq2[j - 1] + optimal_seq2
                i -= 1
                j -= 1
            elif score_matrix[i][j] == score_matrix[i - 1][j - 1] + m:
                optimal_seq1 = seq1[i - 1] + optimal_seq1
                optimal_seq2 = seq2[j - 1] + optimal_seq2
                i -= 1
                j -= 1
            elif score_matrix[i][j] == score_matrix[i - 1][j] + g:
                optimal_seq1 = seq1[i - 1] + optimal_seq1
                optimal_seq2 = "-" + optimal_seq2
                i -= 1
            else:
                optimal_seq1 = "-" + optimal_seq1
                optimal_seq2 = seq2[j - 1] + optimal_seq2
                j -= 1

        alignments.append((optimal_seq1, optimal_seq2))

    return alignments
 

if __name__ == "__main__":
    # Default scoring scheme
    match = 6
    mismatch = -3
    gap = -3

    if len(sys.argv) == 3: 
        match = int(sys.argv[1])
        mismatch = int(sys.argv[2])
        gap = int(sys.argv[3])

    s1 = "" 
    s2 = ""
    s1, s2 = parse_FASTA_file(sys.stdin) 
    alignments = local_align(s1, s2, match, mismatch, gap)

    for i, alignment in enumerate(alignments):
        print(f"Alignment {i+1}: {alignment[0]}")
        print(f"             {alignment[1]}")
