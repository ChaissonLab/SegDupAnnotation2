# SegDupAnnotation2 General Configuration

This workflow requires a `config/sd_analysis.json` configuration file. Please modify it to your needs. It is valudated against `workflow/schemas/config.schema.json` which also contains example parameters.

## Required Parameters

The key files required by this pipeline are an assembly file, a genemodel (for instance a refseq annotation file), and the PacBio reads bam.

| Parameter | Type | Description |
| --- | --- | ------ |
| species | String | The name, code, or identifier for the species or individual being processed. |
| reads_bam | Array of Strings | List of PacBio read file paths in bam format. |
| asm | String | File path to assembly of this species/individual's genome in fasta format. |
| genemodel | String | Fasta filepath of chosen gene model. Assumes fasta headers have no space characters. |
| temp | String | File path to an existing directory for temporary files. |

## Optional Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | ------ |
| read_type | String | N/A | Note PacBio read technology type (CLR vs CCS). Currently used for metadata purposed only. |
| sex_chr | Array of Strings | N/A | List of sex chromosome names in given assembly. |
| min_copy_identity | Number | 0.90 | Minimum gene copy identity to keep when the gene copy is compared to the original copy. |
| min_hit_length | Int | 5000 | Minimum hit length to keep in bases. |
| max_length_margin | Number | 0.10 | Keep gene copies with length within \<max_length_margin\> of the original gene's length. |
| min_depth | Number | 0.05 | Minimum mean copy depth to keep as percentage of mean assembly depth. |
| flag_filter_single_exon_genes | Boolean | true | When true, keeps only genes with multiple exons. |
| flag_assume_clear_and_unique_gene_codes | Boolean | true | When false, assumes gene model fasta headers are in default RefSeq format, and thus renames all headers based on gene symbol in parenthesis at end of header line. |
| flag_filt_uncharacterized_genes | Boolean | true | When true, filters out genes in gene model with gene names beginning with `LOC`. |
| flag_allow_overlapping_genes | Boolean | true | When false, group overlapping genes using network based approach. |
| flag_filtered | Boolean | false | Retain and calculate depth for all genes even those that don't meet filter minimums. |
| grid_xlarge | String | N/A | Generic SLURM/cluster paramtetrs for running certain rules that are computationally intense. Recommended: 64 cores & 130 Gb RAM. |
| grid_large | String | N/A | Generic SLURM/cluster paramtetrs for running certain rules that are computationally intense. Recommended: 16 cores & 48 Gb RAM. |
| grid_medium | String | N/A | Generic SLURM/cluster paramtetrs for running certain rules that are computationally mild. Recommended: 4 cores & 6 Gb RAM. |
| grid_small | String | N/A | Generic SLURM/cluster paramtetrs for running certain rules that are not computationally intense. Recommended: 1 cores & 1 Gb RAM. |
