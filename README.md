# tetranucleotide clustering
This is a simple Perl script for clustering based on tetranucleotide frequency.

It will get the normalised count for the tetranucleotides then calculate the Euclidean distance between each sequences and cluster them based on an unsupervised clustering method, agglomerative clustering.

I have included an optimisation method, the Silhouette coefficient, but it will take many times the time and resources to finish.

## How to use
There are three params:
- i: the input file (single lined fasta)
- a: automatically get the cutoff threshold of the dendrogram and produce the clusters (set 1 to use and 0 to turn off)
- p: use palindromes or not (set 1 to use the palindromes and 0 to use all available tetranucleotides)
- c: if a is turned off then manually set the cutoff here (must be set to a value if a is off)

Example:

perl tetranucleotides.pl -i input.fasta -a 1 ==> cluster based on palindromic tetranucleotides and automated threshold optimisation

perl tetranucleotides.pl -i input.fasta -a 0 -p 0 -c 1.3 ==> cluster based on all tetranucleotides and manually set threshold
