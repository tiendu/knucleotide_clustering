# k-nucleotide clustering
This is a simple Perl script for clustering based on k-nucleotide frequency.

It will get the normalised count for the k-nucleotides then calculate the Euclidean distance between each sequences and cluster them based on an unsupervised clustering method, hierarchical agglomerative clustering.

I have included an optimisation method, the Silhouette coefficient, but it will take many times the time and resources to finish.

## How to use
There are few params:
-i or --input: the input file (single lined fasta)
-k or --knucleotide: as in k-nucleotide e.g., if k = 4 then we are using tetranucleotides for clustering
-a or --auto: automatically get the cutoff threshold of the dendrogram and produce the clusters (set 1 to use and 0 to turn off)
-p or palindromes: use palindromes or not (set 1 to use the palindromes and 0 to use all available tetranucleotides)
-c or cutoff: if a is turned off then manually set the cutoff here (must be set to a value if a is off)

Example:

perl knucleotides.pl -i input fasta -k 4 ==> cluster based on only palindromic tetranucleotides and automated threshold optimisation

perl knucleotides.pl -i input.fasta -k 4 -p 0 ==> cluster based on all tetranucleotides and automated threshold optimisation

perl knucleotides.pl -i input.fasta -k 4 -a 0 -p 0 -c 1.3 ==> cluster based on all tetranucleotides and manually set threshold
