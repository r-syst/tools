# tools
Contains some tools useful in barcoding and metabarcoding   

## disseq

C program computing pairwise distance between a set of sequences from the score of Smith-Waterman or Needleman-Wunsch alignment.

## declic

Python program  for analysing patterns of correlations between molecular based and morphological based taxonomy in reference databases. The inputs are
* the fasta file
* the character file
* the distance file   
as provided in directories of referece databases.   

Main available analysis are:
* Multidimensional Scaling to visualize in low dimenson the structre of a point cloud with distances as close as possible to the genetc distanes in distance file
* Visualize connected components of graphs built from threshholds in distances (barcoding gap).

**Warning** declic is in alpha state: there still are some bugs and it is available for tests and bug reports. 

**contact** alain.franc@inria.fr








