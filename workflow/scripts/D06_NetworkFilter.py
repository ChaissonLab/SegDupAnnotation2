#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 05/24/23

# Purpose: Identify consensus isoforms from isoform clusters using lpa community detection on hits whose exons overlap.
# Input: Bed12 filepath
# Output: First 8 cols of Bed12 input file with a final column noting 0 if filtered out or 1 if hit remains in final set.

### NOTE: At this time, input cannot be piped in, since input file is accessed multiple times.

import sys
import argparse
import numpy as np
import networkx as nx
import pandas as pd

# NETWORK TUNING PARAMS:
MIN_OVERLAP = 50 # remove edges in exon overlap graph if exon overlap is less than this parameter in bases

# Parse Input
print("---- Running D06_NetworkFilter.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Filter out extra isoforms.")
parser.add_argument("bed_filepath", help="Input bed12 file. (Including comma separated list of exon lengths, then a comma separated list of exon start positions.)")
parser.add_argument("communities_filepath", help="Output tab separated file of consensus isoform. (The columns are: consensus isoform, comma separated list of isoforms in the consensus isoform's community.)")
args = parser.parse_args()

# Determine array dimensions
bedFile=open(args.bed_filepath)

GENE_BED_I=3
CHROM_BED_I=0
START_BED_I=1
END_BED_I=2
EXON_SIZES_BED_I=10
EXON_STARTS_BED_I=11

numTotalCopies=0
numTotalExons=0
geneStrLen=0
chromStrLen=0

for line in bedFile:
    bedLine=line.split()
    numTotalCopies+=1
    thisGeneStrLen=len(bedLine[GENE_BED_I])
    if geneStrLen < thisGeneStrLen:
        geneStrLen = thisGeneStrLen
    thisChromStrLen=len(bedLine[CHROM_BED_I])
    if chromStrLen < thisChromStrLen:
        chromStrLen = thisChromStrLen
    numExons=len(bedLine[EXON_STARTS_BED_I].split(","))
    numTotalExons+=numExons

bedFile.close()


# Create ndarrays
isoforms = np.empty([numTotalCopies],
                    dtype={'names': ('chrom', 'start', 'end', 'gene', 'score', 'strand',
                                     'first_exon_start', 'last_exon_end',
                                     'copy_id','picked','exon_sum'),
                          'formats':('U'+str(chromStrLen),'u8','u8','U'+str(geneStrLen),'u2','U1',
                                     'u8','u8',
                                     'u8','u1','u8')})

GENE_ISOFORMS_I=3
COPY_ID_ISOFORMS_I=8
PICKED_ISOFORMS_I=9
EXON_SUM_ISOFORMS_I=10

exons = np.empty([numTotalExons],
                 dtype={'names': ('copy_id', 'chrom', 'start', 'end','sum'),
                        'formats': ('u8','U'+str(chromStrLen),'u8','u8','u8')})

COPY_ID_EXONS_I=0
CHROM_EXONS_I=1
START_EXONS_I=2
END_EXONS_I=3
SUM_EXONS_I=4


# Fill ndarrays
bedFile=open(args.bed_filepath)

copyNum=0
exonNum=0
for line in bedFile:
    bedLine=line.split()
    exonSizes=[int(num) for num in bedLine[EXON_SIZES_BED_I].split(",")]
    exonSum=sum(exonSizes)
    isoforms[copyNum]=(bedLine[CHROM_BED_I],int(bedLine[START_BED_I]),int(bedLine[END_BED_I]),bedLine[GENE_BED_I],
                       int(bedLine[4]),bedLine[5],int(bedLine[6]),
                       int(bedLine[7]),copyNum,0,exonSum)
    exonStarts=bedLine[EXON_SIZES_BED_I].split(",")
    for i in range(len(exonSizes)):
        exonStart=int(bedLine[START_BED_I])+int(exonStarts[i]) # gene start plus exon start (exon start is relative to first position of gene, not chromosome)
        exonEnd=exonStart+int(exonSizes[i])
        exons[exonNum]=(copyNum,bedLine[CHROM_BED_I],exonStart,exonEnd,exonSum) # copy_id, chrom, start, end, sum
        exonNum+=1
    copyNum+=1

bedFile.close()

exons.sort(order=('chrom','start','end','copy_id'))


# Create Graph
print("-- Creating Graph", file=sys.stderr)
# (nodes = index of each gene, edges = exon overlap, weight overlap = num bases of overlap between the exons of two genes)
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
        avg=((isoforms[u][EXON_SUM_ISOFORMS_I] + isoforms[v][EXON_SUM_ISOFORMS_I]) / 2)
        overlapOverAvg = o['overlap'] / avg
        sqOverlapOverAvg = overlapOverAvg * overlapOverAvg # square of average of sum of exons of nodes connected by edge
        overlapGraph.add_edge(u,v,overlap=o['overlap'],
                              overlapOverAvgLen=overlapOverAvg,
                              overlapOverAvgLen_min1=1-overlapOverAvg,
                              sqOverlapOverAvgLen=sqOverlapOverAvg)

del edgeList

print("-  "+str(len(overlapGraph.nodes))+" gene isoforms detected (nodes).", file=sys.stderr)
print("-  "+str(len(overlapGraph.edges))+" isoform to isoform exon overlaps (edges).", file=sys.stderr)

# Community Detection
print("-- Detecting Communities", file=sys.stderr)
subgraphs=[]
communities_lpa=[]
for subgraph_nodes in sorted(nx.connected_components(overlapGraph),key=len,reverse=True):
    sub=overlapGraph.subgraph(subgraph_nodes).copy()
    subgraphs.append(sub)
    coms=nx.community.label_propagation_communities(sub)
    for c in coms:
        communities_lpa.append(c)
print("-  "+str(len(communities_lpa))+" gene communities detected.", file=sys.stderr)

# communities_lpa=sorted(communities_lpa, key=len, reverse=True) # Useful for testing only - TODO Delete Me

coms_out = open(args.communities_filepath,"w")
print("\t".join(["consensus_isoform","all_isoforms_in_community"]),file=coms_out)
for c in communities_lpa:
    # Identify Consensus Node
    big_n = next(iter(c))
    big_oa = 0
    sub_c = overlapGraph.subgraph(c).copy()
    iso_names = []
    for n in c:
        sum_oa = 0
        for u,v in sub_c.edges(n): # for each edge of node n
            sum_oa += sub_c.get_edge_data(u,v)['overlapOverAvgLen']
        if sum_oa > big_oa:
            big_oa = sum_oa
            big_n = n
        iso_names.append(isoforms[n][GENE_ISOFORMS_I])
    # Record Consensus Node
    picked=big_n
    isoforms[picked][PICKED_ISOFORMS_I]=1
    print("\t".join([isoforms[picked][GENE_ISOFORMS_I],
                     ",".join(iso_names)]),
          file=coms_out)
coms_out.close()

# Output Data
print("-- Printing Output", file=sys.stderr)
df=pd.DataFrame(isoforms)
df=df.drop('copy_id',axis=1) # delete copy_id col
df=df.drop('exon_sum',axis=1) # delete exon_sum col

df.to_csv(sys.stdout,sep="\t",encoding='utf-8',header=False,index=False)
