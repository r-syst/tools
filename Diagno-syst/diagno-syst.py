#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__  = "Alain Franc; alain.franc@inra.fr"
__coll__    = "Jean-Marc Frigerio"
__version__ = "0.0.1"
__date__    = "December 13, 2016"

"""
------------------------------------------------------------------------

How it is called

------------------------------------------------------------------------

usage	diagno-syst.py -s <sample> -c <yamlfile> 

@ s		string	name of the sample
@ c		string	name of the configuration file (yaml)

What it does:   builds an inventory by mapping an environmental sample
                on a reference database
                
recommended pratice:
- put diagno-syst.py in bin directory
- start the programm from directory /src of the project
	where the yaml config file should be

------------------------------------------------------------------------

What is read from configuration file

------------------------------------------------------------------------

-----------
Input files
-----------
@ fasfile   : fasta file (NGS) of the environmental sample
@ charfile  : the character file of the references
@ disfile   : the distance file between all queries and all references, sparse format

----------
Parameters
----------
@ charname  : taxonomic level on which to label queries (a header of the character file)
@ gap_max   : upper limit for barcoding gap
@ tab       : if True, the array of inventories, one gap per column, is written
@ inventory : if True, the inventory derived from <Tab> is written (largest nb of reads per column)
@ sort      : if True, taxa in inventory sorted in decreasing order of nb of reads
@ stacks    : if True, the stacks per taxon will be written in a file per taxon in <stacks_taxon.fas>

------------
Output files
------------
@ rep       : directory where to write the output files
              the names of the output files are automatically generated
              and written or not according to the values of
                @ tab           the array of inventories, one gap per column,
                @ inventory     the inventory derived from <Tab>
                @ stacks        the stacks query reads per taxon, a file per taxon


------------------------------------------------------------------------

Notes

------------------------------------------------------------------------

-> How it proceeds:

0 - loads the fasta file
1 - loads the character file, and builds
2 - loads the distance file, built by mpi_disseq, sparse format, asDistance file is in a sparse format
3 - annotates and builds the inventory
4 - Organizes the inventories as a table
5 - Reduce phase: a single inventory

-> With more details ...

0 - loads the fasta file
    dedat library is not used, as it is too slow
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

last revision:  Dec 8, 2016
"""


import numpy as np
import matplotlib.pyplot as plt
import argparse
import os, sys
import yaml

# ----------------------------------------------------------------------

def print_header():
    """
    """
    print("\n--------------------------------------------------------------------")
    print("Supervised clustering")
    print("builds an inventory by mapping an environmental sample")
    print("on a reference database")
    print("")
    print("Version: 0.0.1c")
    print("Date: Dec 9, 2016 - revised March, 22, 2017 ; May 17, 2018")
    print("Author: Alain Franc and Jean-Marc Frigerio")
    print("Maintainer: Alain Franc")
    print("Contact: alain.franc@inra.fr")
    print('Please cite: Frigerio & al. - 2016 - diagno-syst: a tool for')
    print('accurate inventories in metabarcoding - ArXiv:1611.09410')
    print("download @ http://arxiv.org/abs/1611.09410")
    print("--------------------------------------------------------------------\n")

# ----------------------------------------------------------------------

def informative(q, gap, nei, ref_id, ref_taxa):
    """
    what it does: labels a query from the knowledge of its neighborhood
                    method of informative read (the nieghborhood must be of a single color)

    Arguments
    @ q         :
    @ gap       : maximum distance accepted for being in same taxon (barcoding gap)
    @ nei       : neighborhood of the query
    @ ref_id    : seq_id of the reference
    @ ref_taxa  : names of references in character table

    Result

    What it does
    - reads the sparse distance file
    - extracts the three columns: query id's, reference id's, and distance between them
    - then makes a loop over all barcoding gaps
    - for ,each gap
    - looks at all neighbors
    - if no neighbor: the query is labelled "unclassified"
    - if one neighbor: the query is labelled with the taxon of tyhe unique neighbor

    af, last revision 25/10/2016
    """
    l_nei   = len(nei)      # number of neighbors
    if l_nei == 0:
        tax_annot   = "unclassified"   # note: <unknown> may lead to unsolvable problems with some databases afterwards
    else:
        if l_nei==1:        #
            i           = nei[0]
            item        = ref[i]
            if ref_id.count(item) > 0:
                k           = ref_id.index(item)
                tax_annot   = ref_taxa[k]				
            else:
                print("!!warning: item", item, "is not in reference database")  
                tax_annot = "unclassified"          

        else:
            nei_taxa    = l_nei*[""]
            nei_ref     = [ref[i] for i in nei]
            for i, item in enumerate(nei_ref):
                if ref_id.count(item) > 0:
                    k           	= ref_id.index(item)
                    nei_taxa[i] 	= ref_taxa[k]
                else:
                    print("!!warning: item", item, "is not in reference database")
            nei_taxa_div	= list(set(nei_taxa))
            if len(nei_taxa_div)==1:
                tax_annot = nei_taxa[0]
            else:
                tax_annot = "ambiguous"
    return tax_annot

# ----------------------------------------------------------------------

def annotation(queries, dis, ref_id, ref_taxa, gap):
    """
    @ queries   : query identifiers (first column of Dis)
    @ ref_id    : reference identifiers (second column of Dis)
    @ dis       : distance between both
    @ ref_taxa  : name of references in character table
    @ gap       : barcoding gap

    Result
    query_annot : dictionnary   : key is query id, value is annotation
    """
    #query_annot ={}             # dictionary with annotations
    query       = queries[0]
    i           = 0
    j           = i
    n           = len(queries)
    while i < n:
        query   = queries[i]    # query: the query id to be annotated
        j       = i
        # selects all rows with a same query, = from i to j (included)
        while queries[j]==query:
            if j == n-1:
                break
            else:
                j = j+1
        #
        query_rows  = range(i,j+1)                  # rows with a given query
        dis_rows    = [dis[i] for i in query_rows]  # distances with all reference at distace < 30
        which       = [i for i, item in enumerate(dis_rows) if item <= gap] # those at distance <= gap
        nei         = [query_rows[i] for i in which]                        # which row index
        tax_annot   = informative(query, gap, nei, ref_id, ref_taxa)        # query annotation
        if tax_annot not in ['ambiguous','unclassified']:
            query_annot[query] = tax_annot                                  # added to dictonary query_annot
        i           = j+1
    #
    return query_annot

# ----------------------------------------------------------------------

def get_inventory(query_annot, op_print=True):
    """
    @ query_annot   : a dictionary with annotation per query, from annotations

    Result
    inventory       : a dictionary, with key: annotations, values: the number of
                                    queries with this annotation
    """
    taxa_list   = [query_annot[read] for read in query_annot.keys()]
    taxa_names  = list(set(taxa_list))
    inventory   = {}
    for name in taxa_names:
        inventory[name] = taxa_list.count(name)
    #
    if op_print:
        for name in sorted(inventory.keys()):
            print(name + "-> " + str(inventory[name]))
    return inventory
    
def get_options():
    """Parse options"""
    opt_parser = argparse.ArgumentParser(description=__doc__)
    opt_parser.add_argument('-c', '--config', required=True, help='yaml configuraton file')
    opt_parser.add_argument('-s', '--sample', required=True, help='sample under study within the project')
    #
    return opt_parser.parse_args()    


# ======================================================================
#
#                   main()
#
# ======================================================================

print_header()

#------------- Parameters ----------------------------------------------
# imported from params.py

opts 	= get_options()
config	= opts.config
sample 	= opts.sample

param   = yaml.load(open(config))

project     = param['project'] ; print("project directory is", project)
rep         = param['rep'] + project
if param['fasfile'] !="None":
	fasfile     = rep + param['fasfile']  + sample + '.fas'
else:
	fasfile=None	
disfile     = rep + param['disrep']  + sample + param['dis_suffix']
charfile    = rep + param['charfile']
charname    = param['charname']
gap_max     = param['gap_max']      # maximal gap (estimated: 3% of the length)
sort        = param['sort']      # taxa in inventory sorted in decreasing order
tab         = param['tab']    # table of inventories for all gaps
inventory   = param['inventory']    # writes the inventory on the disk
stacks      = param['stacks']

#raw_input()
# ----------------------------------------------------------------------


if fasfile:
    # 0 - loading fasta file
    #
    #   note: dedat library is not used, as it is too slow
    #         the fasta file must by in following format:
    #           - even line (0,2,4, ...): >seq_id
    #           - odd line  (1,3,5, ...): the sequence, without carriage return

    print("loading fasta file ...")
    print("fasfile = ", fasfile)
    fasta       = np.loadtxt(fasfile, dtype=object)
    n           = len(fasta)
    query_seqid = [fasta[i][1:] for i in range(0, n, 2)]
    query_word  = [fasta[i] for i in range(1, n, 2)]
    query_len   = [len(word) for word in query_word]
    #plt.hist(query_len, 50, normed=1, facecolor='g', alpha=0.75)
    #plt.show()
    print("fasta file loaded")



# 1 - loading character file
#       ref_taxa    : list of taxa of <charname> in character file
#       ref_id      : corresponding seq_id

print("loading character file ...")
Char            = np.loadtxt(charfile, delimiter='\t', comments='#', dtype=str)
print("Character file loaded from ")
print(charfile)
char_headers    = list(Char[0,:])
k               = char_headers.index(charname)
ref_taxa        = Char[1:,k]
ref_id          = list(Char[1:,0])


# 2 - loading distance file
#
# Distance file is in a sparse format
# Dis has three columns
# they are extracted here
#   - first     : queries   : query id
#   - second    : ref       : reference id
#   - third     : dis       : distance between both

print("loading distance file; this may take a while ...")
Dis         = np.loadtxt(disfile, delimiter="\t", dtype=str)
n           = Dis.shape[0]                      # number of rows in DIS
queries     = Dis[:,0]                          # query id's
query_set   = list(set(queries))                # list of query id's (uniq)
ref         = Dis[:,1]                          # references which matched per query
dis         = np.array(Dis[:,2],dtype=float)    # distance of match
print("distance file loaded from ", disfile)

# 3 - annotation and inventory
#
# pl_inv    : a list of inventories, one per gap
#             each inventory is a dictionary
#             with    the key being the taxons in reference
#                     the values being the number of reads with this annotation, fot this gap

print("annotating ...")
pl_inv      = []
# Initialization of query_annot
#   dictionary with annotations
#   keys:   queries seq_id
#   values: annotation
query_annot = {}
for  item in queries:
    query_annot[item]   = "unclassified"
#
for gap in range(1+gap_max):
    """
    makes a loop for each gap between 0 and gap_max (included)
    """
    print("gap = " + str(gap) + " ..." + chr(13),end='')
    query_annot = annotation(queries, dis, ref_id, ref_taxa, gap)   # annotates all queries
    inv         = get_inventory(query_annot, op_print=False)        # inventory as a dictionary
    pl_inv.append(inv)
print()


# 4 - Organizing the inventories as a table
taxa    = []        # a list which will have all taxa present in at least one inventory
for a in range(len(pl_inv)):
    item        = pl_inv[a]
    new_taxa    = list(item.keys())
    taxa        = taxa + new_taxa
taxa    = list(set(taxa))
n_taxa  = len(taxa)
taxa.sort()
# Res is a table with taxa as rows, gaps as columns, and Res[i,j] is the number of reads of taxon <i> for gap <j>
Res = np.zeros(shape=(n_taxa, 1+gap_max), dtype=int)
for gap in range(1+gap_max):
    dic = pl_inv[gap]
    for item in dic.keys():
        k               = taxa.index(item)
        Res[k, gap]     = dic[item]


# if a resfile, writes Res as this file
if tab:
    resfile = rep + "inventories/" + sample  + "_inventory_per_gap.txt"
    os.makedirs(rep + 'inventories', exist_ok=True)
    with open(resfile, "w") as out_handle:
        gaps    = ["gap_"+str(i) for i in range(1+gap_max)]
        headers = ["Taxon"] + gaps
        item    = "\t".join(headers)
        out_handle.write(item + '\n')
        for i in range(n_taxa):
            vals    = [str(item) for item in Res[i,:]]
            head    = [taxa[i]]
            row     = head + vals
            item = "\t".join(row)
            out_handle.write(item + '\n')
        out_handle.close()
    print("\nInventories per gap written in", resfile, "\n")

# 5 - Reduce phase: a single inventory

# the number of reads for the inventory is the maximum number over all gaps
nb_reads    = Res.max(axis=1)

# if sort == True, sorts the taxa in decreasing number of number of reads
if sort:
    ind         = np.argsort(nb_reads)
    ind         = [ind[i] for i in np.arange(start = len(ind)-1, stop=-1, step=-1)]
    nb_reads    = [nb_reads[i] for i in ind]
    taxa        = [taxa[i] for i in ind]

# prints the inventory on screen
print("\nInventory with taxa in decreasing order ...")
for a, item in enumerate(taxa):
    print(item + " -> " + str(nb_reads[a]))


# if inventory, writes the inventory in file inventory
if inventory:
	os.makedirs('../inventories',exist_ok=True)
	inventoryfile   = "../inventories/" + sample + ".inv"
	with open(inventoryfile, "w") as out_handle:
		headers = ["Taxon", "Nb"]
		item    = "\t".join(headers)
		out_handle.write(item + '\n')
		print('TAXA', taxa)
		for a, item in enumerate(taxa):
			if item not in ['ambiguous','unclassified']:
				row     = [item, str(nb_reads[a])]
				item    = "\t".join(row)
				out_handle.write(item + '\n')
		out_handle.close()
	print("\nInventory written in", inventoryfile, "\n")

if stacks:
    print("\n-------------------------------------")
    print("writing the fasta files per taxon ...")
    print("-------------------------------------\n")
    for item in taxa:
        if item not in ['ambiguous','unclassified']:
            print(item)
            seqfile = rep + "stacks/" + item + ".fas"
            with open(seqfile, "w") as out_handle:
                which = [query_id for query_id in query_annot.keys() if query_annot[query_id] == item]
                print(str(len(which)), "reads")
                for q_id in which:
                    if query_seqid.count(q_id)==1:
                        k       = query_seqid.index(q_id)
                        word    = query_word[k]
                        wtw     = ">" + q_id
                        out_handle.write(wtw + '\n')
                        out_handle.write(word + '\n')
                    else:
                        print("query " + q_id + " is not in fasta file ...")
                out_handle.close()
    #

# ======================================================================
#
#           That's all folks!
#
# ======================================================================













