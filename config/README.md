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
| flag_filt_single_exon_genes | Boolean | true | When true, keeps only genes with multiple exons. |
| flag_assume_clear_and_unique_gene_codes | Boolean | true | When false, assumes gene model fasta headers are in default RefSeq format, and thus renames all headers based on gene symbol in parenthesis at end of header line. |
| flag_filt_uncharacterized_genes | Boolean | true | When true, filters out genes in gene model with gene names beginning with `LOC`. |
| flag_allow_overlapping_genes | Boolean | true | When false, group overlapping genes using network based approach. |
| flag_filtered | Boolean | false | Retain and calculate depth for all genes even those that don't meet filter minimums. |
| cluster_mem_mb_baby | Int | 1000 | The memory in MB a cluster node or cpu must provide for a computationally simple job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_small | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally simple job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_medium | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally mild job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_large | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally intense job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_xlarge | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally intense job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_baby | Int | 1 | The number of cpus per task for a computationally simple rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_small | Int | 3 | The number of cpus per task for a computationally mild rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_medium | Int | 3 | The number of cpus per task for a computationally intense rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_large | Int | 3 | The number of cpus per task for a computationally intense rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_runtime_short | Int | 240 | The walltime in minutes allocated for rules expected to take a relatively short amount of time (like 4 hrs). This parameter is only used if called in the cluster_exec parameter or by snakemake's --slurm command line paramter. |
| cluster_runtime_long | Int | 240 | The walltime in minutes allocated for rules expected to take a relatively long amount of time (like 24 hrs). This parameter is only used if called in the cluster_exec parameter or by snakemake's --slurm command line paramter. |
| override_mem | Int | -1 | Override the memory available in MB otherwise defined by the cluster_mem_mb_\<size\> parameters in MB. If set to -1, the cluster_mem_mb_<X> paramters will not be overwritten. |
| override_num_cores | Int | -1 | Override the number of allocated cores otherwise defined by the cluster_cpus_per_task_\<size\> parameters. If set to -1, the cluster_cpus_per_task_\<size\> paramters will not be overwritten. |
