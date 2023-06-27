#!/usr/bin/env python

# Input:
#   Arg1: bam of resolved originals (see rule D02_FilterResolvedOriginals_FiltPercentAligned in snakefile).  
# Purpose: Filter out input bam entries with <= 50% match (between assembly and reads).

import pysam
import sys

resOrig = pysam.AlignmentFile(sys.argv[1])

sys.stdout.write(str(resOrig.header))

for aln in resOrig.fetch():
    astat=aln.get_cigar_stats()
    nHardClip=astat[0][5]
    nAligned = astat[0][0] + astat[0][7] + astat[0][8] # = alignment match + sequence match + sequence mismatch
    
    if (aln.query_sequence is None):
        continue
    qLen = len(aln.query_sequence) + nHardClip

    if nAligned / qLen >= 0.5:
        sys.stdout.write(aln.to_string() + "\n")

resOrig.close()