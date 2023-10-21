#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 05/24/23

# Purpose: Identify consensus isoforms from isoform clusters using lpa community detection on hits whose exons overlap.
# Input: 
# Output: Same file as input with final column noting 0 if filtered out or 1 if hit remains in final set.

### NOTE: At this time, input cannot be piped in, since input file is accessed multiple times.

import sys
import argparse
import numpy as np
import networkx as nx
import pandas as pd

# NETWORK TUNING PARAMS:
MIN_OVERLAP = 50 # remove edges in exon overlap graph if exon overlap is less than this parameter in bases

# Parse Input
parser = argparse.ArgumentParser(description="Filter out all but one isoform for each gene.")
parser.add_argument("pafxe_filepath", help="Input pafxe file. pafxe is first 12 columns from paf format followed by 8 columns defined in script E02_CalcPafIdentity.py. Then followed by comma separated list of exon lengths, then a comma separated list of exon start positions.")
args = parser.parse_args()

# Determine array dimensions to use
pafxeFile=open(args.pafxe_filepath)

GENE_PAFXE_I=0
CHROM_PAFXE_I=5
START_PAFXE_I=7
END_PAFXE_I=8
EXON_SIZES_PAFXE_I=20
EXON_STARTS_PAFXE_I=21

numTotalCopies=0
numTotalExons=0
geneStrLen=0
chromStrLen=0

for line in pafxeFile:
    pafxeLine=line.split()
    numTotalCopies+=1
    thisGeneStrLen=len(pafxeLine[GENE_PAFXE_I])
    if geneStrLen < thisGeneStrLen:
        geneStrLen = thisGeneStrLen
    thisChromStrLen=len(pafxeLine[CHROM_PAFXE_I])
    if chromStrLen < thisChromStrLen:
        chromStrLen = thisChromStrLen
    numExons=len(pafxeLine[EXON_STARTS_PAFXE_I].split(","))
    numTotalExons+=numExons

pafxeFile.close()

# Create ndarrays
pafx = np.empty([numTotalCopies],
                dtype={'names': ('gene', 'q_seq_len', 'q_start', 'q_end', 'strand', 'chrom', 't_seq_len',
                                'start', 'end', 'num_matching_bases_in_mapping', 'num_bases_in_mapping',
                                'mapping_qual', 'aln_type', 'num_chain_minimizers', 'chaining_score',
                                'chaining_score_of_secondary', 'seq_divergence',
                                'length_of_query_rgn_with_repetitive_seeds', 'percent_identity',
                                'percent_accuracy','copy_id','picked','exon_sum'),
                        'formats':('U'+str(geneStrLen),'u8','u8','u8','U1','U'+str(chromStrLen),'i',
                                'u8','u8','u8','u8',
                                'u1','U1','U10','U10',
                                'U10','U10',
                                'U10','d',
                                'd','u8','u1','u8')})

P_IDENTITY_PAFX_I=18
P_ACCURACY_PAFX_I=19
COPY_ID_PAFX_I=20
PICKED_PAFX_I=21
EXON_SUM_PAFX_I=22

exons = np.empty([numTotalExons],
                 dtype={'names': ('copy_id', 'chrom', 'start', 'end','sum'),
                        'formats': ('u8','U'+str(chromStrLen),'u8','u8','u8')})

COPY_ID_EXONS_I=0
CHROM_EXONS_I=1
START_EXONS_I=2
END_EXONS_I=3
SUM_EXONS_I=4


# Fill ndarrays
pafxeFile=open(args.pafxe_filepath)

copyNum=0
exonNum=0
for line in pafxeFile:
    pafxeLine=line.split()
    exonSizes=[int(num) for num in pafxeLine[20].split(",")]
    exonSum=sum(exonSizes)
    pafx[copyNum]=(pafxeLine[0],int(pafxeLine[1]),int(pafxeLine[2]),int(pafxeLine[3]),pafxeLine[4],pafxeLine[5],int(pafxeLine[6]),
                   int(pafxeLine[7]),int(pafxeLine[8]),int(pafxeLine[9]),int(pafxeLine[10]),
                   int(pafxeLine[11]),pafxeLine[12],pafxeLine[13],pafxeLine[14],
                   pafxeLine[15],pafxeLine[16],
                   pafxeLine[17],float(pafxeLine[18]),
                   float(pafxeLine[19]),copyNum,0,exonSum)
    exonStarts=pafxeLine[21].split(",")
    for i in range(len(exonSizes)):
        exonStart=int(pafxeLine[7])+int(exonStarts[i]) # gene start plus exon start (exon start is relative to first position of gene, not chromosome)
        exonEnd=exonStart+int(exonSizes[i])
        exons[exonNum]=(copyNum,pafxeLine[5],exonStart,exonEnd,exonSum) # copy_id, chrom, start, end, sum
        exonNum+=1
    copyNum+=1

pafxeFile.close()

exons.sort(order=('chrom','start','end','copy_id'))


# Create Graph (nodes = index of each gene, edges = exon overlap, weight overlap = num bases of overlap between the exons of two genes)
overlapGraph = nx.Graph()
# Add Nodes
for geneNum in range(numTotalCopies):
    overlapGraph.add_node(geneNum)

# Helper functions for adding edges.
# Function overlap: returns true if two genes overlap and false 
# otherwise. Assumes start position of first gene is less than
# start position of second gene.
# c1=chrom_gene1, e1=end_gene1, c2=chrom_gene2, s2=start_gene2
def overlap (c1,e1,c2,s2):
    return (c1==c2 and s2<=e1)

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
for ex in exons:
    j=i+1
    while j<exons.size and overlap(ex[CHROM_EXONS_I],ex[END_EXONS_I],exons[j][CHROM_EXONS_I],exons[j][START_EXONS_I]):
        if overlapGraph.has_edge(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I]):
            overlapSum=overlapGraph.get_edge_data(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I])["overlap"] + overlapDist(ex[END_EXONS_I],exons[j][START_EXONS_I],exons[j][END_EXONS_I])
            overlapGraph.add_edge(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I],overlap=overlapSum)
        else:
            overlapGraph.add_edge(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I],overlap=overlapDist(ex[END_EXONS_I],exons[j][START_EXONS_I],exons[j][END_EXONS_I]))
        j+=1
    i+=1

# Calculating Heuristic Weights for each edge and pruning short overlaps
edgeList=list(overlapGraph.edges(data=True))
for u,v,o in edgeList:
    if o['overlap'] < MIN_OVERLAP:
        overlapGraph.remove_edge(u,v)
    else:
    # average of sum of exons of nodes connected by edge
        avg=((pafx[u][EXON_SUM_PAFX_I] + pafx[v][EXON_SUM_PAFX_I]) / 2)
        overlapOverAvg = o['overlap'] / avg
        sqOverlapOverAvg = overlapOverAvg * overlapOverAvg # square of average of sum of exons of nodes connected by edge
        overlapGraph.add_edge(u,v,overlap=o['overlap'],
                              overlapOverAvgLen=overlapOverAvg,
                              overlapOverAvgLen_min1=1-overlapOverAvg,
                              sqOverlapOverAvgLen=sqOverlapOverAvg)

del edgeList

# Community Detection
subgraphs=[]
communities_lpa=[]
for subgraph_nodes in sorted(nx.connected_components(overlapGraph),key=len,reverse=True):
    sub=overlapGraph.subgraph(subgraph_nodes).copy()
    subgraphs.append(sub)
    coms=nx.community.label_propagation_communities(sub)
    for c in coms:
        communities_lpa.append(c)

# communities_lpa=sorted(communities_lpa, key=len, reverse=True) # TODO probs unnecessary

for c in communities_lpa:
    # Identify Consensus Node
    big_n = next(iter(c))
    big_oa = 0
    sub_c=overlapGraph.subgraph(c).copy()
    for n in c:
        sum_oa=0
        for u,v in sub_c.edges(n): # for each edge of node n
            sum_oa += sub_c.get_edge_data(u,v)['overlapOverAvgLen']
        if sum_oa > big_oa:
            big_oa = sum_oa
            big_n = n
    # Record Consensus Node
    picked=big_n
    pafx[picked][PICKED_PAFX_I]=1

# Print pafxe # TODO No longer necessary
# pafx.sort(order=('copy_id'))

df=pd.DataFrame(pafx)
df=df.drop('copy_id',axis=1) # delete copy_id col
df=df.drop('exon_sum',axis=1) # delte exon_sum col

df.to_csv(sys.stdout,sep="\t",encoding='utf-8',header=False,index=False)
