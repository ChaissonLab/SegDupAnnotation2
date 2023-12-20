#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 11/01/23

# Purpose: Pick consensus label

import sys
import argparse
import numpy as np
import networkx as nx
import pandas as pd

print("---- Running E10_PickConsensusLabel.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Determine Consensus Labels for Isoform Communities.")
parser.add_argument("communities_tsv_filepath", help="Input tsv file. (Only 2nd col used. 2nd col should have comma separated list of gene names.)")
parser.add_argument("communities_consol_tsv_filepath", help="tsv file created with col1 having represenative label and col2 has a comma separated list of associated labels in its community.")
parser.add_argument("bed_in_filepath", help="Bed filepath with gene name in col4 to be replaced by representative family name.")
parser.add_argument("bed_out_filepath", help="Identical to input bed filepath but gene labels are switched with representative labels.")
args = parser.parse_args()

# Create Label Graph
print("-- Creating Graph", file=sys.stderr)
labelGraph = nx.Graph()

# Add Nodes and Edges
comsFile=open(args.communities_tsv_filepath)
for line in comsFile:
    comsLine=line.split()
    labels=comsLine[1].split(sep=",")
    for i in range(len(labels)):
        labelGraph.add_node(labels[i])
        for j in range(i+1,len(labels)):
            labelGraph.add_edge(labels[i],labels[j]) # adds nodes if not in graph and edge if not in graph
comsFile.close()

# Arbitrarily pick representative node from graph and create communities output file
print("-- Picking Representative Node", file=sys.stderr)
label_dict={} # Create Dictionary # every label has a key. It's value is the picked/representative label of whose family it belongs to.
coms_out = open(args.communities_consol_tsv_filepath,'w')
#coms_out = open("results/E10_rep_communities.tsv",'w')
for subgraph_nodes in sorted(nx.connected_components(labelGraph),key=len,reverse=True):
    # subgraph_nodes is a set
    picked=next(iter(subgraph_nodes))
    label_dict[picked] = picked
    print(picked+'\t'+','.join(list(subgraph_nodes)),file=coms_out)
    for n in subgraph_nodes:
        label_dict[n] = picked
coms_out.close()

# Regenerate Bed input with representative gene label names
print("-- Updating bed file", file=sys.stderr)
bedIn=open(args.bed_in_filepath)
bedOut=open(args.bed_out_filepath,'w')
#bedIn=open("results/E09_resolved_copies_geneLabelsNotGrouped.bed")
#bedOut=open("results/E10_resolved_copies.bed",'w')
for line in bedIn:
    vals=line.split()
    vals[3]=label_dict[vals[3]]
    print('\t'.join(vals),file=bedOut)
bedIn.close()
bedOut.close()



