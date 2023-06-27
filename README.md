# SegDupAnnotation2

## Overview

SegDupAnnotation2 is a snakemake workflow designed to count gene duplications given PacBio reads, genome assembly, and a gene model.  
Successor to [SegDupAnnotation](https://github.com/ChaissonLab/SegDupAnnotation).

## System Requirements

TODO

## Installation

Clone github repository. Create configuration file per config/README.md specifications. Then run `snakemake -c 1 -j 250 --use-conda` or  `snakemake -c 1 -j 250 --cluster "{params.grid_opts}" --use-conda -k`.

## Salient Output File Specifications

results/G01\_dups.bed
| Column | Description |
| --- | ------ |
| #chr | Gene copy's position in assembly. |
| start | |
| end | |
| gene | Gene name. |
| orig\_chr | Position of gene copy's original copy in assembly. |
| orig\_start |  |
| orig\_end | |
| strand | Strand on which gene copy is on: 0 for 'Original', '1' for reverse. |
| p\_identity | Similarity to 'Original' gene calculated as: #matches/(#matches+#mismatches+#insertion\_events+#deletion\_events) |
| p\_accuracy | Similarity to 'Original' gene calculated as : #matches/(#matches+#mismatches+#insertions+#deletions) |
| identity | Notes whether copy is the 'Original' or a resolved 'Copy'. |
| depth | Mean gene depth over mean assembly depth calculated using 100 bp bins. |
| depth\_stdev | Standard deviation of gene depth over mean assembly depth calculated using 100 bp bins. |
| copy\_num | Rounded depth value. |
| depth\_by\_vcf | Copy number as determined by hmm's vcf output. |

results/G03\_per\_gene\_counts.bed
| Column | Description |
| --- | ------ |
| gene | Gene name. |
| depthNormalizedToAsm | Mean of depth across all gene copies divided by mean assembly depth. |
| depth | Sum of depth across all gene copies. |
| measuredCov | depth rounded to nearest whole number. |
| copyCount | Number of total gene copies - resolved plus collapsed. |
| resolvedCount | Number of resolved gene copies. |

results/G04\_summary\_stats.tsv
Contains tab delimited summary statistics.

## Tips

If rerunning pipeline using the same input assembly and reads, retain and copy over files A01\_assembly.fasta and A04\_assembly.bam as aligning reads to assembly is the longest step.
