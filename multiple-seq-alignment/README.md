# gibbs_sampler.py

#### 1. To run program from command line:
```
python3 gibbs_sampler.py < input.fasta > output.txt
```
where [input.fasta] is any FASTA-formatted file containing multiple sequences. 

#### 2. Program overview:
gibbs_sampler.py
- Parses the streamed input FASTA file into a list of multiple string sequences.
- Runs iterations of the Gibbs sampling method through the class SimpleGibbsSampler until a modified convergence criterion is satisfied.
- Returns an MSA of the current best motif instance in each sequence.
