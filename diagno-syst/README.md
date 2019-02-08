# Diagno-syst

Python program for supervized clustering of environmental samples

*version 0.1*    
*Authors: Alain Franc & Jean-Marc Frigerio*    
*Maintainer: Alain Franc*    
*contact: alain.franc@inra.fr*    


## Usage   


``diagno-syst.py -s <sample> -c <yamlfile>`` 


with  

* ``-s``	string	name of the sample	  
* ``-c``	string	name of the configuration file (yaml)   


### What it does

Builds an inventory by mapping an environmental sample on a reference database, by supervized clustering   

### Quick install

``diagon-syst.py`` is a stand alone python program which requires some python dependencies:
* numpy
* matplotlib
* yaml   

The installation is given here for Linux. Once you have installed python 3.6 or higher
* ``pip install numpy``
* ``pip install matplotlib``
* ``pip install pyyaml``    


There are other possibilities, like installing Anaconda, which includes numpy and matplotlib, but not yaml.    


Running ``diagno-syst.py`` requires that the following files are available:
* a distance file (distance between queries and references, computed beforehand)
* a parameter file, specifying the input files, output files, and selected options(all options and files in one single file)
* a reference database, containing a taxonomic annotation of the reference sequences used for distances.   


To install ``diagno-syst.py``, and the files used in the example (a distance file, a parameter file, a reference database)
* type ``git clone https://github.com/r-syst/tools.git``
This will install all tools (not only ``diagno-syst``), with one directory per tool. This creates a subdirectory ``tools`` where you are, and a subdirectory ``diagno-syst`` 
within ``tools``. Then, go into this subdirectory by
* ``cd tools/diagno-syst``   

and you are in the directory where the program and the data for running the examples are located. To get all the data, you have to uncompress the distance file. For this
* type ``bunzip2 SN1-55_027_dissw.txt.bz2``
* type ``chmod +x diagno-syst.py`` (just to be sure ...)   

You are now ready to run the example.



### Quick start on an example

After installaton as above, type in  terminal in ``[...]/tools/diagno-syst$``
* ``./diagno-syst.py -c params_rbcL.yaml -s SN1-55_027``   


This creates a subdirectory ``inventories`` in directory ``diagno-syst``, and (according to the parameters in configuraton file) writes two inventories:
* one oer gap
* one synthetic over all selected gaps



## What is read from configuration file

### Input files

* ``fasfile``   : fasta file (NGS) of the environmental sample
* ``charfile``  : the character file of the references
* ``disfile``   : the distance file between all queries and all references, sparse format

### Parameters

* ``charname``  : taxonomic level on which to label queries (a header of the character file)
* ``gap_max``   : upper limit for barcoding gap
* ``tab``       : if True, the array of inventories, one gap per column, is written
* ``inventory`` : if True, the inventory derived from <Tab> is written (largest nb of reads per column)
* ``sort``      : if True, taxa in inventory sorted in decreasing order of nb of reads
* ``stacks``    : if True, the stacks per taxon will be written in a file per taxon in <stacks_taxon.fas>


## Notes


### How it proceeds:

0 - loads the fasta file   
1 - loads the character file   
2 - loads the distance file, built by mpi_disseq, sparse format
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

2 - loads the distance file, built by mpi_disseq, sparse format, as an array <Dis> with three columns   
    - first     : queries   : query id   
    - second    : ref       : reference id   
    - third     : dis       : distance between both   

3 - annotates and builds the inventory    
  
    pl_inv  = []    : a list of inventories, one per gap   
                      each inventory is a dictionary   
                      with    keys    : the taxons in reference   
                              values  : the number of reads with this annotation, for this gap   

    query_annot = {}    dictionary with annotations   
                        keys    : queries (see [2])   
                        values  : a taxon if any, or "unclassified"    

4 - Organizes the inventories as a table   

    Res is a table with   
    * row         : taxa    
    * colums      : gaps   
    * Res[i,j]    : the number of reads of taxon <i> for gap <j>   
    
    if Tab == True   
        Res is written as <rep_fulldiagno.txt>   

5 - Reduce phase: a single inventory   
    the number of reads for the inventory is the maximum number over all gaps   
