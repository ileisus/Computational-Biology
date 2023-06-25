""" 
GlobalSequenceAlignment

Summary
-------
Performs global DNA sequence alignment using a similarity scoring scheme.
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


class GlobalAlign:
    """ Implementation of the Needleman-Wunsch Algorithm 

    Example Usage
    -------------
    # Align 2 sequences using a defined scoring scheme 
    align = GlobalSeqAlignment(seq1, seq2, match, mismatch, gap)

    # Get optimal global sequence alignments
    x_prime, y_prime = align.get_optimal_alignment()

    Attributes
    ----------
    s1: sequence 1
    s2: sequence 2
    M: match score
    m: mismatch penalty
    g: gap penalty
    """

    def __init__(self, s1, s2, M, m, g):
        self.s1 = s1   
        self.s2 = s2
        self.M = M   
        self.m = m   
        self.g = g 

    def match_score(self, s1, s2):
        if s1 == s2:
            return self.M
        elif s1 == '-' or s2 == '-':
            return self.g
        else:
            return self.m

    def optimal_align_helper(self, s1, s2): 
        """ 
        """
        m = len(s1)  
        n = len(s2)

        score = np.zeros((m+1, n+1)).astype(int).tolist()

        # Calculate score table
        for i in range(0, m+1):
            score[i][0] = self.g * i
        for j in range(0, n+1):
            score[0][j] = self.g * j
        for i in range(1, m+1):
            for j in range(1, n+1):
                match = score[i-1][j-1] + self.match_score(s1[j-1], s2[i-1]) 
                delete = score[i-1][j] + self.g
                insert = score[i][j-1] + self.g
                score[i][j] = max(match, delete, insert)

        # Traceback
        align1 = ""
        align2 = ""
        i = m
        j = n
        
        while i > 0 and j > 0:
            score_curr = score[i][j]
            score_diag = score[i-1][j-1]
            score_up = score[i][j-1]
            score_left = score[i-1][j]

            if score_curr == score_diag + self.match_score(s1[j-1], s2[i-1]):
                align1 += s1[j-1]
                align2 += s2[i-1]
                i -= 1
                j -= 1
            elif score_curr == score_up + self.g:
                align1 += s1[j-1]
                align2 += '-'
                j -= 1
            elif score_curr == score_left + self.g:
                align1 += '-'
                align2 += s2[i-1]
                i -= 1

        while j > 0:
            align1 += s1[j-1]
            align2 += '-'
            j -= 1
        while i > 0:
            align1 += '-'
            align2 += s2[i-1]
            i -= 1

        # Compensate for reverse traceback order
        align1 = align1[::-1]
        align2 = align2[::-1]

        return(align1, align2)

    def get_optimal_alignment(self):
        """
        Returns
        -------
            Two strings containing optimal sequence alignments
        """
        return(self.optimal_align_helper(self.s1, self.s2))


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
    align = GlobalAlign(s1, s2, match, mismatch, gap)
    s1_prime, s2_prime = align.get_optimal_alignment()
    print(s1_prime + "\n" + s2_prime)

