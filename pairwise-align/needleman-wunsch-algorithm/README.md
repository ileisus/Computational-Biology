#needleman-wunsch

# 1. Summary 

The Needleman-Wunsch algorithm is often used to perform pairwise, global sequence alignment. 

needleman-wunsch.py reads in two DNA sequences stored in a FASTA file, computes global alignments, and returns the optimal alignment.
Traceback preference = diagonal > upward > leftward traceback.
Scoring scheme is defined by the user upon execution.


# 2. Run/execute 
- Run executable with command
      python align.py [M_award] [m_penalty] [g_penalty] < [someFASTA.py]
            
  where
      M_award = user-defined match score,
      m_penalty = user-defined mismatch penalty,
      g_penalty = user-defined gap penalty,
      someFASTA.py = user-specified sequence input file.


