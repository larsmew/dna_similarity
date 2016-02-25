Application to compare DNA sequencing datasets.
==============
Features:
- Support for identification of Single-Nucleotide Polymorphisms (SNPs)

The method uses two steps:

  1. Locality-sensitive hashing to find candidate pairs of similar reads.
  2. Candidate pairs are used as input to construct groups (or blocks) of overlapping reads in which substitution mutations are identified.

Application options:

      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -n <FILENAME>, --normal=<FILENAME>
                            Fasta file containing normal reads.
      -d <FILENAME>, --diseased=<FILENAME>
                            Fasta file containing diseased (mutated) reads.
      -i <FILENAME>, --candidate_pairs=<FILENAME>
                            set <FILENAME> as input file containing
                            candidate pairs to import.
      -e <FILENAME>, --export_candidate_pairs=<FILENAME>
                            set <FILENAME> as output file to export
                            candidate pairs to. Dump as either txt, json
                            or pickle file. Just put the right extension
                            and everything else is automated.
      -l <FILENAME>, --log_file=<FILENAME>
                            set <FILENAME> as log file.
      -k <VALUE>, --k_size=<VALUE>
                            set <VALUE> as size for k-mers.
                            Default: 12
      -b <Value>, --bands=<Value>
                            set <VALUE> as the number of bands for LSH.
                            Default: 25
      -r <Value>, --rows=<Value>
                            set <VALUE> as the number of rows for LSH.
                            Default: 4
      -s <VALUE>, --seed=<VALUE>
                            set <VALUE> as seed for hash functions.
      -S <VALUE>, --supporting_reads=<VALUE>
                            Reads required to support mutation.
                            Default: 3
      -o <VALUE>, --mismatch_overlap=<VALUE>
                            Number of allowed mismatches in overlap.
                            Default: 1
      -g <VALUE>, --mismatch_group=<VALUE>
                            Number of allowed mismatches in group.
                            Default: 2


Example usage:
  
      >>> python compareReads.py -n normal_reads.fasta -d diseased_reads.fasta -k 12 -b 2 -r 5 -l log.txt -S 3 -o 1 -g 2
      
      Computing table of shingles positions...      
      Finished computation of shingles table in 6.61666666667e-06 minutes
      Number of k-mers: 47
      Memory usage (in mb): 49.8515625
      
      Minhashing...
      Finished minhashing in 0.000328633333333 minutes
      Memory usage (in mb): 50.703125
      Running LSH and finding similar pairs...
      Number of buckets in band 1: 4
      Number of candidate pairs in band 1: 42
      Number of unique candidate pairs in band 1: 24
      Minhashing...
      Finished minhashing in 0.000119683333333 minutes
      Memory usage (in mb): 50.8359375
      Running LSH and finding similar pairs...
      Number of buckets in band 2: 3
      Number of candidate pairs in band 2: 74
      Number of unique candidate pairs in band 2: 40
      
      Number of unique candidate pairs 42
      Finished LSH in 0.000569666666667 minutes
      Memory usage (in mb): 50.83984375 
      
      Collecting reads from file fasta files
      Finished reading in 2.81666666667e-06 minutes
      Found 5 sequences in fasta files
      Memory usage (in mb): 50.83984375
      
      read_R: 3
      Consensus:
       GGAAAACCCCTTTTGGGGCCCCTTTTAAAAGGAAAAGAAACCCCCCCCTTTTTTTTGGGGGGGGAAAACCCCTTTTGGGGCCCCTTTTAAAAGG
                                                                   C                                 
       GGAAAACCCCTTTTGGGGCCCCTTTTAAAAGGAAAAGAAACCCCCCCCTTTTTTTTGGGGGG*
                                       AAAAGAAACCCCCCCCTTTTTTTTGGGGGGGGAAAACCCCTTTTGGGGCCCCTTTTAAAAGG
                                      GAAAAGAAACCCCCCCCTTTTTTTTGGGGGGGGAAAACCCCTTTTGGGGCCCCTTTTAAAAG
        GAAAACCCCTTTTGGGGCCCCTTTTAAAAGGAAAAGAAACCCCCCCCTTTTTTTTGGGGGGG
      
                                       AAAAGAAACCCCCCCCTTTTTTTTGGGGCGGGAAAACCCCTTTTGGGGCCCCTTTTAAAAGG
                                      GAAAAGAAACCCCCCCCTTTTTTTTGGGGCGGGAAAACCCCTTTTGGGGCCCCTTTTAAAAG
       GGAAAACCCCTTTTGGGGCCCCTTTTAAAAGGAAAAGAAACCCCCCCCTTTTTTTTGGGGCG
        GAAAACCCCTTTTGGGGCCCCTTTTAAAAGGAAAAGAAACCCCCCCCTTTTTTTTGGGGCGG
      Left parts: [0, 6, 10, 16]
      Number of left parts: 4
      Right parts: [5, 13, 15]
      Number of right parts: 3
      mismatches: [60]
