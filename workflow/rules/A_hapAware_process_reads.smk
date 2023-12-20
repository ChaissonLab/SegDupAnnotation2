# HAPLOTYPE AWARE MODE ACTIVE
rule A01H_linkAsm:
    input:
        hap1=config["asm"],
        hap2=config["hap2"]
    output:
        hap1link="results/A01H_hap1_asm.fasta",
        hap2link="results/A01H_hap2_asm.fasta"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/A01H_linkAsm.log"
    shell:"""
        echo "##### A01H_linkAsm" > {log}
        ln -s {input.hap1} {params.workflowDir}/../{output.hap1link} 2>> {log}
        ln -s {input.hap2} {params.workflowDir}/../{output.hap2link} 2>> {log}
        """

rule A02H_faiIndexAsm:
    input:
        hap="results/A01H_hap{hapNum}_asm.fasta"
    output:
        hapfai="results/A01H_hap{hapNum}_asm.fasta.fai",
        haplnk="results/A02H_hap{hapNum}_asm.fai"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A02H_faiIndexAsm.{hapNum}.log"
    benchmark: "benchmark/A02H_faiIndexAsm.{hapNum}.tsv"
    shell:"""
        echo "##### A02H_faiIndexAsm - hap{wildcards.hapNum}" > {log}
        echo "### Create Samtools fai Index"
        samtools faidx {input.hap} --fai-idx {output.hapfai} 2>> {log}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{output.hapfai} {params.workflowDir}/../{output.haplnk}
    """

rule A03H_mmiIndexAsm:
    input:
        hap="results/A01H_hap{hapNum}_asm.fasta"
    output:
        mmi="results/A03H_hap{hapNum}_asm.mmi"
    params:
        workflowDir=workflow.basedir,
        read_type=config["read_type"]
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A03H_mmiIndexAsm.hap{hapNum}.log"
    benchmark: "benchmark/A03H_mmiIndexAsm.hap{hapNum}.tsv"
    shell:"""
        echo "##### A03H_mmiIndexAsm - hap{wildcards.hapNum}" > {log}
        echo "### Determine MM2 Parameters" >> {log}
        if [ {params.read_type} = "CCS" ]
        then
            mm2_mode="map-hifi"
        else
            mm2_mode="map-pb"
        fi
        echo "Minimap2 -x preset: $mm2_mode" >> {log}

        echo "### Create Minimap2 mmi Index" >> {log}
        minimap2 -x $mm2_mode -d {output.mmi} {input.hap} -t {resources.cpus_per_task} 2>> {log}
    """

rule A04H_CombineHaps:
    input:
        hap1="results/A01H_hap1_asm.fasta",
        hap2="results/A01H_hap2_asm.fasta"
    output:
        asm="results/A04H_asm.fasta",
        fai="results/A04H_asm.fasta.fai"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A04H_CombineHaps.log"
    benchmark: "benchmark/A04H_CombineHaps.tsv"
    shell:"""
        echo "##### A04H_CombineHaps" > {log}
        echo "### Concat Haps" >> {log}
        cat {input.hap1} {input.hap2} > {output.asm}

        echo "### Index Asm" >> {log}
        samtools faidx {output.asm} --fai-idx {output.fai} 2>> {log}
    """

rule A05H_alignReads:
    input:
        bam=lambda wildcards: bamFiles[wildcards.base],
        mmi=ancient("results/A03H_hap{hapNum}_asm.mmi")
    output:
        fastq=temp("results/A05H_hap{hapNum}_aligned/A05H_{base}.fastq"),
        aligned=temp("results/A05H_hap{hapNum}_aligned/A05H_{base}.bam")
    priority: 10
    params:
        read_type=config["read_type"]
    resources:
        tmpdir=tmpDir,
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A05H_alignReads.hap{hapNum}.{base}.log"
    benchmark: "benchmark/A05H_alignReads.hap{hapNum}.{base}.tsv"
    shell:"""
    {{
        echo "##### A05H_alignReads - hap{wildcards.hapNum} - bam: {wildcards.base}" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)

        echo "### Determine MM2 Parameters" >> {log}
        if [ {params.read_type} = "CCS" ]
        then
            mm2_mode="map-hifi"
        else
            mm2_mode="map-pb"
        fi
        echo "Minimap2 -x preset: $mm2_mode" >> {log}

        echo "### Make Tmp Dir" >> {log}
        tmp_collate_path=`mktemp -d -p {resources.tmpdir} A05H.collate.XXXX.tmp`
        tmp_mm2_path=`mktemp -d -p {resources.tmpdir} A05H.mm2.XXXX.tmp`
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} A05H.sort.XXXX.tmp`

        echo "### Align Reads" >> {log}
        samtools view -h -F 2304 -u -@ "$numAdditionalThreads" {input.bam} | \
            samtools collate -O -u -@ "$numAdditionalThreads" - "$tmp_collate_path" | \
            samtools fastq -@ "$numAdditionalThreads" - > {output.fastq}
        minimap2 {input.mmi} {output.fastq} -a -x "$mm2_mode" -t {resources.cpus_per_task} --split-prefix "$tmp_mm2_path" | \
            samtools sort -T "$tmp_sort_path"/reads -m "$mem_per_cpu"M -@ "$numAdditionalThreads" -o {output.aligned}

        echo "### Delete Tmp Dirs and their Contents" >> {log}
        rm -rf "$tmp_collate_path" "$tmp_mm2_path" "$tmp_sort_path"
    }} 2>> {log}
    """ # TODO in minimap use -K5M ? Also probably don't need collate step on raw reads

rule A06H_mergeReads:
    input:
        aln=expand("results/A05H_hap{{hapNum}}_aligned/A05H_{base}.bam", base=bamFiles.keys())
    output:
        mrg=protected("results/A06H_hap{hapNum}_asm_indludesAutozygousRgns.bam")
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A06H_mergeReads.hap{hapNum}.log"
    benchmark: "benchmark/A06H_mergeReads.hap{hapNum}.tsv"
    shell:"""
    {{
        echo "##### A06H_mergeReads - hap{wildcards.hapNum}" > {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)
        samtools merge {output.mrg} {input.aln} -@"$numAdditionalThreads"
    }} 2>> {log}
    """

rule A07H_baiIndexBam:
    input:
        bam="results/A06H_hap{hapNum}_asm_indludesAutozygousRgns.bam"
    output:
        bai="results/A06H_hap{hapNum}_asm_indludesAutozygousRgns.bam.bai",
        lnk="results/A07H_hap{hapNum}_asm_indludesAutozygousRgns.bai"
    conda: "../envs/sda2.main.yml"
    log: "logs/A07H_baiIndexBam.hap{hapNum}.log"
    benchmark: "benchmark/A07H_baiIndexBam.hap{hapNum}.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
    {{
        echo "##### A07H_baiIndexBam - hap{wildcards.hapNum}" > {log}
        echo "### Index Bam" >> {log}
        samtools index {input.bam} -@{resources.cpus_per_task}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{input.bam}.bai {params.workflowDir}/../{output.lnk}
    }} 2>> {log}
    """

rule A08H_partition_reads:
    input:
        #mrg=protected("results/A06H_hap{hapNum}_asm_indludesAutozygousRgns.bam")
        hap1bam="results/A06H_hap1_asm_indludesAutozygousRgns.bam",
        hap2bam="results/A06H_hap2_asm_indludesAutozygousRgns.bam"
    output:
        hap1_reads=temp("results/A08H_hap1_reads.txt"),
        hap2_reads=temp("results/A08H_hap2_reads.txt"),
        hap1_read_names=temp("results/A08H_hap1_read_names.txt"),
        hap2_read_names=temp("results/A08H_hap2_read_names.txt"),
        hapX_read_names=temp("results/A08H_autozygous_read_names.txt"),
        tsv="results/A08H_reads_partitioned.tsv",
        hap1_final_bam="results/A08H_hap1_reads_only.bam",
        hap2_final_bam="results/A08H_hap2_reads_only.bam",
        hapX_final_bam="results/A08H_autozygous_reads_only.bam"
    conda: "../envs/sda2.main.yml"
    log: "logs/A08H_partition_reads.log"
    benchmark: "benchmark/A08H_partition_reads.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
    {{
        echo "##### A08H_partition_reads" > {log}
        tmp_sort_path_hap1=`mktemp -d -p {resources.tmpdir} A08H_hap1.sort.XXXX.tmp`
        tmp_sort_path_hap2=`mktemp -d -p {resources.tmpdir} A08H_hap2.sort.XXXX.tmp`

        echo "### Sort Reads by Name" >> {log}
        samtools view {input.hap1bam} | \
            cut -f1 | \
            sort -u --parallel={resources.cpus_per_task} -T "$tmp_sort_path_hap1" > {output.hap1_reads}
        
        samtools view {input.hap2bam} | \
            cut -f1 | \
            sort -u --parallel={resources.cpus_per_task} -T "$tmp_sort_path_hap2" > {output.hap2_reads}
        
        echo "### Partition Read Names" >> {log}
        comm {output.hap1_reads} {output.hap2_reads} > {output.tsv}
        cat {output.tsv} | cut -f1 | sed '/^[[:space:]]*$/d' > {output.hap1_read_names}
        cat {output.tsv} | cut -f2 | sed '/^[[:space:]]*$/d' > {output.hap2_read_names}
        cat {output.tsv} | cut -f3 | sed '/^[[:space:]]*$/d' > {output.hapX_read_names}

        echo "### Create Read Partitioned bams" >> {log}
        samtools view -b {input.hap1bam} -N {output.hap1_read_names} > {output.hap1_final_bam}
        samtools view -b {input.hap2bam} -N {output.hap2_read_names} > {output.hap2_final_bam}
        samtools view -b {input.hap1bam} -N {output.hapX_read_names} > {output.hapX_final_bam}

        echo "### Delete Tmp Dirs and their Contents" >> {log}
        rm -rf "$tmp_sort_path_hap1" "$tmp_sort_path_hap2"
    }} 2>> {log}
    """

rule A_temp:
    output:
        t3=touch("results/A01_assembly.fasta"),
        t4=touch("results/A01_assembly.fasta.fai"),
        t1=touch("results/A04_assembly.bam"),
        t2=touch("results/A05_assembly.bai"),
        hapX_final_bam="results/A06_autozygous_reads_only.bam"
    localrule: True
    shell: """ 
    echo "TODO DELETE THIS RULE"
    """


# rule A07_align_partitioned_reads_to_hap1_union_hap2:
#     input:
#         hap1bam="results/A04_assembly_hap1.bam",
#         hap2bam="results/A04_assembly_hap2.bam"
#     output:
#         hap1_reads=temp("A06_hap1_reads.txt"),
#         hap1_reads=temp("A06_hap1_reads.txt"),
#         tsv="results/A06_reads_partitioned.tsv",
#         hap1_final_bam="results/A06_hap1_reads_only.bam",
#         hap2_final_bam="results/A06_hap2_reads_only.bam",
#         hapX_final_bam="results/A06_autozygous_reads_only.bam"
#     conda: "../envs/sda2.main.yml"
#     log: "logs/A07_align_partitioned_reads_to_hap1_union_hap2.log"
#     benchmark: "benchmark/A07_align_partitioned_reads_to_hap1_union_hap2.tsv"
#     params:
#         workflowDir=workflow.basedir
#     resources:
#         mem_mb=cluster_mem_mb_small,
#         cpus_per_task=cluster_cpus_per_task_small,
#         runtime=config["cluster_runtime_long"]
#     shell:"""
#     {{
#         echo "##### A07_align_partitioned_reads_to_hap1_union_hap2 - Haplotype Aware" > {log}
#         samtools view {input.hap1bam} -N A06_hap1_read_names.txt > {output.hap1_final_bam}
#         samtools view {input.hap2bam} -N A06_hap2_read_names.txt > {output.hap2_final_bam}
#         samtools view {input.hap1bam} -N A06_hapX_read_names.txt > {output.hapX_final_bam}

#     }} 2>> {log}
#     """