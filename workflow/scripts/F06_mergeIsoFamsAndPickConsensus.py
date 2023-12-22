#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 12/21/23

# Purpose: Take grouped isoforms table and merge groups that contain the same isoforms.
# Input: Grouped isoforms table.
# Output: Grouped isoforms table with overlapping groups merged.

import sys
import argparse
import pathlib
import networkx as nx

# Parse Input
print("---- Running F06_mergeIsoFamsAndPickConsensus.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Take grouped isoforms table and merge groups that contain the same isoforms.")
parser.add_argument("tsv_filepath", type=pathlib.Path, help="Input tab separated file of grouped isoforms. (The first column is a comma separated list of isoform codes without its location information. The second column has a comma separated list of isoform codes with each copy's location information.) (Only the first column is needed for this script.)")
args = parser.parse_args()

# Create Graph
# (nodes = gene name, edges = conections given in input tsv_filepath file)
gene_fam_graph = nx.Graph()

# Parse Input File
tsvFile=open(args.tsv_filepath)
for line in tsvFile:
    tsvLine=line.split()
    genes=tsvLine[0].split(',')
    for g in genes:
        gene_fam_graph.add_edge(genes[0],g)
tsvFile.close()

gene_fam_names = nx.connected_components(gene_fam_graph)

for fam_nodes in gene_fam_names:
    subgraph = gene_fam_graph.subgraph(fam_nodes).copy()
    print(",".join(list(subgraph.nodes())), file=sys.stdout)
