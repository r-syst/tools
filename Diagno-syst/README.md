# Diagno-syst

Python program for supervized clustering of environmental samples

*version 0.1*    
*Authors: Alain Franc & Jean-Marc Frigerio*    
*Maintainer: Alain Franc*    
*contact: alain.franc@inra.fr*    


## Usage   


``diagno-syst.py -s <sample> -c <yamlfile>`` 


with```

* -s		string	name of the sample	  
* -c		string	name of the configuration file (yaml)   


### What it does

Builds an inventory by mapping an environmental sample on a reference database, by supervized clustering   

### recommended pratice

* put diagno-syst.py in bin directory
* start the programm from directory (which can be named) /src of the project where the yaml config file should be


## What is read from configuration file

### Input files

* ``fasfile``   : fasta file (NGS) of the environmental sample
* ``charfile``  : the character file of the references
* ``disfile``   : the distance file between all queries and all references, sparse format

### Parameters

* charname  : taxonomic level on which to label queries (a header of the character file)
* gap_max   : upper limit for barcoding gap
* tab       : if True, the array of inventories, one gap per column, is written
* inventory : if True, the inventory derived from <Tab> is written (largest nb of reads per column)
* sort      : if True, taxa in inventory sorted in decreasing order of nb of reads
* stacks    : if True, the stacks per taxon will be written in a file per taxon in <stacks_taxon.fas>

### Output files

* rep       : directory where to write the output files   
              the names of the output files are automatically generated   
              and written or not according to the values of   
                @ tab           the array of inventories, one gap per column,   
                @ inventory     the inventory derived from <Tab>    
                @ stacks        the stacks query reads per taxon, a file per taxon   

## Notes


### How it proceeds:

0 - loads the fasta file   
1 - loads the character file, and builds   
2 - loads the distance file, built by mpi_disseq, sparse format, asDistance file is in a sparse format    
3 - annotates and builds the inventory   
4 - Organizes the inventories as a table   
5 - Reduce phase: a single inventory   

### With more details ...    

0 - loads the fasta file   
    the fasta file must by in following format:   
        - even line (0,2,4, ...): >seq_id   
        - odd line  (1,3,5, ...): the sequence, without carriage return   

1 - loads the character file, and builds   
    ref_taxa    the taxa of column <charname> in character file   
    ref_id      the corresponding reference seq_id   

2 - loads the distance file, built by mpi_disseq, sparse format, asDistance file is in a sparse format    
    an array <Dis> with three columns   
    - first     : queries   : query id   
    - second    : ref       : reference id   
    - third     : dis       : distance between both   

3 - annotates and builds the inventory   
    pl_inv  : a list of inventories, one per gap   
                each inventory is a dictionary   
                with    keys    : the taxons in reference   
                        values  : the number of reads with this annotation, for this gap   

    query_annot = {}    dictionary with annotations   
                        keys    : queries (see [2])   
                        values  : a taxon if any, or "unclassified"    

4 - Organizes the inventories as a table   
    Res is a table with   
        row         : taxa    
        colums      : gaps   
        Res[i,j]    : the number of reads of taxon <i> for gap <j>   
    if Tab == True   
        Res is written as <rep_fulldiagno.txt>   

5 - Reduce phase: a single inventory   
    the number of reads for the inventory is the maximum number over all gaps   
