"""
Performs multiple sequence alignment and motif discovery using 
Gibbs sample methods. 
"""

import sys
import numpy as np
import math

def parse_FASTA_file(input):
    s = "" 
    for line in input:
        if line.startswith(">"): 
            s += "!"
        else:
            stripped = line.rstrip("\n") 
            s += stripped

    list = s.split("!")
    list.pop(0) 

    return(list)


class SimpleGibbsSampler:
    """
    A class used to perform Gibbs Sampling on a set of sequences
    
    Attributes
    ----------
    seqs: list of str 
        a list of formatted strings representing sequences
    motif_len: int
        the fixed motif length
    lang: list of str 
        the alphabet, i.e. list of the nucleotide or protein language symbols
    lang_prob: float
        the background probability of each language symbol 
    seed: int
        a seed value for random numebr generation
    """

    def __init__(self, seqs, motif_len, lang, lang_prob, seed):
        self.seqs = seqs
        self.motif_len = motif_len
        self.lang = lang
        self.lang_prob = lang_prob
        self.random_state = np.random.RandomState(seed)

    def pick_init_positions(self):
        """ Returns a list of random starting positions for each sequence."""
        start_pos_list = np.zeros_like((self.seqs), dtype = int)
        for i in range(len(self.seqs)):
            start_pos_list[i] = self.random_state.randint(0, len(self.seqs))
            
        return(start_pos_list)
        
    def build_pseudo_count_matrix(self, motifs):
        """ Returns a counts matrix with pseudocount 1 
        Args
        ----
        motifs: A list of sequence strings
        """
        count_matrix = np.ones((len(self.lang), self.motif_len))        
        
        for i in range(len(motifs)):              
            for j in range(self.motif_len):   
                for k in range(len(self.lang)): 
                    if motifs[i][j] == self.lang[k]:
                        count_matrix[k][j] += 1  
                                
        return(count_matrix)
        
    def build_pssm(self, count_matrix):
        """ Returns a log-odds position-specific scoring matrix (PSSM).
        The PSSM is a numpy array of shape (length of alphabet, motif length).
        Args
        ----
        count_matrix: a counts matrix, numpy array
        """
        column_sum = len(self.lang) + len(self.seqs) - 1
        freq_matrix = count_matrix/column_sum 
        
        log_func = lambda t: math.log2(t/self.lang_prob) 
        lgfunc = np.vectorize(log_func) 
        self.pssm = lgfunc(freq_matrix)
        
        return(self.pssm)

    def score_seq_windows(self, seq, pssm):
        """ Returns a list of log-odds scores for each window (i.e. motif)
        The PSSM is a numpy array of shape (length of alphabet, motif length).
        Args
        ----
        count_matrix: numpy array
            a counts matrix with pseudoscounts of 1 
        """
        windows = [] 
        for i in range(len(seq) - self.motif_len + 1):
            windows.append(seq[i: (i + self.motif_len)]) 

        # Stores scores for position at matching index.
        scores_list = np.zeros_like(windows, dtype=float)
 
        for index in range(len(windows)):           
            temp_score = 0                       
            string = windows[index]              
            for pos in range(self.motif_len):         
                for i in range(len(self.lang)):         
                    if string[pos] == self.lang[i]:
                        temp_score += pssm[i][pos]
            scores_list[index] = temp_score

        return(scores_list)
        
    def get_msa(self): 
        """ 
        Returns a list of motifs. Includes k motifs for the k sequences.
        """
        msa = [] 
        for i in range(len(self.seqs)):
            site = self.init_pos[i]
            curr = self.seqs[i][site: site + self.motif_len]
            msa.append(curr)

        return(msa)
    
    def run_sampler(self):
        """ 
        Returns an MSA of the current best motif instance in each sequence.
        """
        last_heldout_index = -99999 # Updated every iteration
        itr_tracker = 0
        converged = False

        # 1. Randomly pick a set of k initial positions for k sequences
        self.init_pos = self.pick_init_positions()
        
        # Repeat until {I*} converges for some heldout sequence
        while converged == False:
            itr_tracker += 1
            print("Iteration number: " + str(itr_tracker)) 
            # 2. Randomly pick a heldout sequence 
            curr_heldout_index = self.random_state.randint(0, len(self.seqs))
            # Ensure no immediate repeats for heldout sequence 
            while (curr_heldout_index == last_heldout_index):
                curr_heldout_index = self.random_state.randint(0, len(self.seqs)) 
                if curr_heldout_index != last_heldout_index:
                    break

            last_heldout_index = curr_heldout_index    
            s_star = self.seqs[curr_heldout_index] 

            # Print out s∗ for current iteration and its index
            print("Heldout s*: " + str(self.seqs[curr_heldout_index]) + " at index " + str(curr_heldout_index) + " in seqs.")
            # 3. Build PSSM for (k-1) seqs 
            curr_msa = self.get_msa() # Get full MSA 

            # Print out MSA with k many items 
            print("The current MSA for the " + str(itr_tracker) + "th iteration is:")
            print(curr_msa)
            
            # Remove heldout seq from MSA, result = k-1 
            curr_msa.pop(curr_heldout_index)
            # Use MSA to build counts, frequency and log-odds matrix 
            count_matrix = self.build_pseudo_count_matrix(curr_msa)
            pssm = self.build_pssm(count_matrix)
            print("The PSSM used for this iteration is: ")
            print(pssm)

            # 4. Score windows of s_star using PSSM 
            scores_for_initpos = self.score_seq_windows(s_star, pssm)

            # 5. Update starting position for s* in the I array 
            max_init_pos = np.argmax(scores_for_initpos)

            # If the best scoring pos for s* is equal to its current starting pos
            if self.init_pos[curr_heldout_index] == max_init_pos:
                converged = True
            else: # Otherwise, update starting pos for s* by corresponding index
                self.init_pos[curr_heldout_index] = max_init_pos
            print("\n")
            if converged == True:
                break
            
        return(self.get_msa())


if __name__ == "__main__":

    # Set default parameters
    motif_len = 6
    lang = ['A', 'C', 'G', 'T']
    lang_prob = 0.25
    seed = 20
    
    seqs = [] 
    seqs = parse_FASTA_file(sys.stdin) 
    sgs = SimpleGibbsSampler(seqs, motif_len, lang, lang_prob, seed)
    msa = sgs.run_sampler() 
