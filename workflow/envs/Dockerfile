FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="6206f51b78be89366a80e485e7b52728a92112ca53c59868f211adeebee3e7c7"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/sda2.hmcnc.yml
#   prefix: /conda-envs/f0c1b13e61342ffac5c089dc49cc57c0
#   name: sda2_hmcnc_env
#   # Conda may fail in solving environment if channel_priority is set to strict.
#   # If it fails, check your ~/.condarc config file's channel_priority param.
#   channels:
#     - bioconda
#     - conda-forge
#     - krabbani
#     - nodefaults
#   dependencies:
#     - bioconda::snakemake>=7.25.0 # 7.25.0
#     - boost # 1.82.0
#     - intervaltree # 3.0.2
#     - tabix # 1.11
#     - hmcnc # 1.0
#   # Version nums commented out following dependencies are versions used by authors.
RUN mkdir -p /conda-envs/f0c1b13e61342ffac5c089dc49cc57c0
COPY workflow/envs/sda2.hmcnc.yml /conda-envs/f0c1b13e61342ffac5c089dc49cc57c0/environment.yaml

# Conda environment:
#   source: workflow/envs/sda2.main.yml
#   prefix: /conda-envs/78e0579661129ebd76bf4ba9ac9f8707
#   name: sda2_main_env
#   # Conda may fail in solving environment if channel_priority is set to strict.
#   # If it fails, check your ~/.condarc config file's channel_priority param.
#   channels:
#     - bioconda
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - snakemake>=7.25.0 # 7.25.0
#     - samtools>=1.11 # 1.17
#     - bedtools # 2.30.0
#     - minimap2 # 2.22
#   # Version nums commented out following dependencies are versions used by authors.
RUN mkdir -p /conda-envs/78e0579661129ebd76bf4ba9ac9f8707
COPY workflow/envs/sda2.main.yml /conda-envs/78e0579661129ebd76bf4ba9ac9f8707/environment.yaml

# Conda environment:
#   source: workflow/envs/sda2.python.yml
#   prefix: /conda-envs/a886ae8a033dbbeb9a37b8166c907fe6
#   name: sda2_python_env
#   # Conda may fail in solving environment if channel_priority is set to strict.
#   # If it fails, check your ~/.condarc config file's channel_priority param.
#   channels:
#     - bioconda
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - snakemake>=7.25.0 # 7.25.0
#     - samtools # 1.12
#     - pysam # 0.17.0
#     - numpy # 1.24.2
#     - pandas # 1.5.3
#     - networkx # 3.1
#     - biopython # 1.79
#     - conda-forge::edlib # 1.3.9
#   # Version nums commented out following dependencies are versions used by authors.
RUN mkdir -p /conda-envs/a886ae8a033dbbeb9a37b8166c907fe6
COPY workflow/envs/sda2.python.yml /conda-envs/a886ae8a033dbbeb9a37b8166c907fe6/environment.yaml

# Conda environment:
#   source: workflow/envs/sda2.r.yml
#   prefix: /conda-envs/a7a4d28aadcd4f144d8dfd3a8ecc712b
#   name: sda2_r_env
#   # Conda may fail in solving environment if channel_priority is set to strict.
#   # If it fails, check your ~/.condarc config file's channel_priority param.
#   channels:
#     - bioconda
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - snakemake>=7.25.0 # 7.25.0
#     - r-base # 4.1.1
#     - r-ggplot2 # 3.4.1
#     - r-gridextra # 2.3
#   # Version nums commented out following dependencies are versions used by authors.
RUN mkdir -p /conda-envs/a7a4d28aadcd4f144d8dfd3a8ecc712b
COPY workflow/envs/sda2.r.yml /conda-envs/a7a4d28aadcd4f144d8dfd3a8ecc712b/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/f0c1b13e61342ffac5c089dc49cc57c0 --file /conda-envs/f0c1b13e61342ffac5c089dc49cc57c0/environment.yaml && \
    mamba env create --prefix /conda-envs/78e0579661129ebd76bf4ba9ac9f8707 --file /conda-envs/78e0579661129ebd76bf4ba9ac9f8707/environment.yaml && \
    mamba env create --prefix /conda-envs/a886ae8a033dbbeb9a37b8166c907fe6 --file /conda-envs/a886ae8a033dbbeb9a37b8166c907fe6/environment.yaml && \
    mamba env create --prefix /conda-envs/a7a4d28aadcd4f144d8dfd3a8ecc712b --file /conda-envs/a7a4d28aadcd4f144d8dfd3a8ecc712b/environment.yaml && \
    mamba clean --all -y
