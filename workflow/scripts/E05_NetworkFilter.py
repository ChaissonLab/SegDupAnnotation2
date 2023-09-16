#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 05/24/23

# Purpose: 
# Input: 
# Output: 

import sys
import argparse
import numpy as np
import networkx as nx
import pandas as pd

# Parse Input
parser = argparse.ArgumentParser(description="Filter out all but one isoform for each gene.")
parser.add_argument("pafx_filepath", help="Input pafx file. Pafx is first 12 columns from paf format followed by 8 columns defined in script E02_CalcPafIdentity.py. Then followed by line number of pafx file, then a zero.")
args = parser.parse_args()

# Read Pafx File
# pafx_ndarray = np.loadtxt(args.pafx_filepath,
#                           dtype={'names': ('gene','chrom','start','end','lineNum','picked'),
#                                  'formats':('a60','a20','d','d','d','i1')}, # TODO may need to (substantially) expand char lengths of strings
#                           usecols=(0,5,7,8,21,22))
pafx_ndarray = np.loadtxt(args.pafx_filepath,
                          dtype={'names': ('gene', 'q_seq_len', 'q_start', 'q_end', 'strand', 'chrom', 't_seq_len',
                                           'start', 'end', 'num_matching_bases_in_mapping', 'num_bases_in_mapping',
                                           'mapping_qual', 'aln_type', 'num_chain_minimizers', 'chaining_score',
                                           'chaining_score_of_secondary', 'seq_divergence',
                                           'length_of_query_rgn_with_repetitive_seeds', 'percent_identity',
                                           'percent_accuracy','lineNum','picked'),
                                 'formats':('a60','i','i','i','a1','a20','i',
                                            'i','i','i','i',
                                            'i','a1','a10','a10',
                                            'a10','a10',
                                            'a10','d',
                                            'd','i','i1')}) # TODO may need to (substantially) expand char lengths of strings

#'q_gene', 'q_seq_len', 'q_start', 'q_end', 'strand', 't_chr', 't_seq_len', 't_start', 't_end', 'num_matching_bases_in_mapping', 'num_bases_in_mapping', 'mapping_qual', 'aln_type', 'num_chain_minimizers', 'chaining_score', 'chaining_score_of_secondary', 'seq_divergence', 'length_of_query_rgn_with_repetitive_seeds', 'percent_identity', 'percent_accuracy', 'edlib_cigar'

pafx_ndarray.sort(order=('chrom','start','end','gene'))

# Visualize array # TODO Delete me
#print(pafx_ndarray)

GENE_I=0
CHROM_I=5
START_I=7
END_I=8
LINE_NUM_I=20
PICKED_I=21

# Create Graph (nodes = index of each gene, edges = overlap)
overlapGraph = nx.Graph()
# Add Nodes
i=0
for gene in pafx_ndarray:
    overlapGraph.add_node(i)
    i+=1

# Helper functions for adding edges.
# Function overlap: returns true if two genes overlap and false 
# otherwise. Assumes start position of first gene is less than
# start position of second gene.
# c1=chrom_gene1, e1=end_gene1, c2=chrom_gene2, s2=start_gene2
def overlap (c1,e1,c2,s2):
    return (c1==c2 and s2<e1)

# Function overlapDist: returns magnitude of overlap in bases of
# two genes. Assumes start position of first gene is less than
# start position of second gene and that both genes overlap.
# e1=end_gene1, s2=start_gene2, e2=end_gene2
def overlapDist (e1,s2,e2):
    if e1 < e2:
        return (e1-s2)
    else:
        return (e2-s2)

# Add Edges
i=0
for gene in pafx_ndarray:
    j=i+1
    while j<pafx_ndarray.size and overlap(gene[CHROM_I],gene[END_I],pafx_ndarray[j][CHROM_I],pafx_ndarray[j][START_I]):
        overlapGraph.add_edge(i,j,overlap=overlapDist(gene[END_I],pafx_ndarray[j][START_I],pafx_ndarray[j][END_I]))
        j+=1
    i+=1

# Print Graph # TODO Delete me
# print("nodes:")
# print(list(overlapGraph.nodes))
# print("edges:")
# print(list(overlapGraph.edges))
# print(len(list(overlapGraph.edges)))
# print(overlapGraph.edges[7,8]['overlap'])
# print(overlapGraph.edges[7,9]['overlap'])
# print(overlapGraph.edges[8,9]['overlap'])

# Community Detection
# print ("######")
communities = nx.community.louvain_communities(overlapGraph) #, seed=823) # TODO For testing
# print(communities)
# print(len(communities))

for c in communities:
    cPopped=c.pop()
    pafx_ndarray[cPopped][PICKED_I]=1
    #sys.stdout.write(str(pafx_ndarray[cPopped])+"\n")

# Print Pafx
pafx_ndarray.sort(order=('lineNum'))

#pafx_ndarray.tofile(sys.stdout,sep='\n',format='%s')

df=pd.DataFrame(pafx_ndarray)
df.drop('lineNum',axis=1) # delete lineNum col

df['gene']=df['gene'].str.decode('utf-8')
df['strand']=df['strand'].str.decode('utf-8')
df['chrom']=df['chrom'].str.decode('utf-8')
df['aln_type']=df['aln_type'].str.decode('utf-8')
df['num_chain_minimizers']=df['num_chain_minimizers'].str.decode('utf-8')
df['chaining_score']=df['chaining_score'].str.decode('utf-8')
df['chaining_score_of_secondary']=df['chaining_score_of_secondary'].str.decode('utf-8')
df['seq_divergence']=df['seq_divergence'].str.decode('utf-8')
df['length_of_query_rgn_with_repetitive_seeds']=df['length_of_query_rgn_with_repetitive_seeds'].str.decode('utf-8')

#sys.stdout.write(df.head(n=5).to_string())
#sys.stdout.write(df.to_string())
#df.head(n=5).to_csv(sep="\t",encoding='utf-8')
df.to_csv(sys.stdout,sep="\t",encoding='utf-8',header=False,index=False)

# pafIn=open(args.pafx_filepath)
# i=0
# for line in pafIn:
#     i+=1
#     vals=line.split()
#     sys.stdout.write(line.rstrip('\n')+'\t'+str(pafx_ndarray[int(vals[LINE_NUM_I])-1][PICKED_I])+"\n")
# pafIn.close()

# if int(vals[LINE_NUM_I]) != i:
#     sys.stderr.write("### Error: lines not re-aligning properly (dimensions off)."+"\n")
#     sys.stderr.write("### i="+str(i)+"; lineNum="+str(vals[LINE_NUM_I])+". Should be equivalent.\n") # Probably tighten this error

    # TODO actually, add new col to np array saying if selected or not, then print whole np array (probably resort array too). Or create tmp files. Or don't read file strait into np array, but build it myself.
    # Probably add some sort of distance of edge of sequence to remove before assessing overlap??


















