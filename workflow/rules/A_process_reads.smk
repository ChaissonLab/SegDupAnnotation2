rule A01_linkAsm:
    input:
        asm=config["asm"]
    output:
        asmlink="results/A01_assembly.fasta"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/A01_linkAsm.log"
    shell:"""
        echo "##### A01_linkAsm" > {log}
        ln -s {input.asm} {params.workflowDir}/../{output.asmlink} 2>> {log}
    """

rule A02_indexAsm:
    input:
        asm="results/A01_assembly.fasta"
    output:
        fai="results/A01_assembly.fasta.fai",
        lnk="results/A02_assembly.fai"
    params:
        cluster_exec=config["cluster_exec"],
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A02_indexAsm.log"
    benchmark: "benchmark/A02_indexAsm.tsv" # TODO Remove all benchmark lines
    shell:"""
        echo "##### A02_indexAsm" > {log}
        samtools faidx {input.asm} --fai-idx {output.fai} 2>> {log}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{output.fai} {params.workflowDir}/../{output.lnk}
    """

rule A03_alignReads:
    input:
        bam=lambda wildcards: bamFiles[wildcards.base],
        asm=ancient("results/A01_assembly.fasta")
    output:
        aligned=temp("results/A03_aligned/A03_{base}.bam")
    priority: 10
    params:
        cluster_exec=config["cluster_exec"]
    resources:
        tmpdir=tmpDir,
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A03_alignReads.{base}.log"
    benchmark: "benchmark/A03_alignReads.{base}.tsv"
    shell:"""
    {{
        echo "##### A03_alignReads" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}

        echo "### Make Tmp Dir" >> {log}
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} reads.sort.XXXX.tmp`

        echo "### Align Reads" >> {log}
        samtools view -h -F 2304 {input.bam} | \
            samtools fastq - | \
            minimap2 {input.asm} - -a -t {resources.cpus_per_task} | \
            samtools sort -T "$tmp_sort_path"/reads -m "$mem_per_cpu"MB -@ $(( {resources.cpus_per_task}-1 )) -o {output.aligned}

        echo "### Delete Samtools Sort Tmp Dir and its Contents" >> {log}
        rm -rf "$tmp_sort_path"
    }} 2>> {log}
    """

rule A04_mergeReads:
    input:
        aln=expand("results/A03_aligned/A03_{base}.bam", base=bamFiles.keys())
    output:
        mrg=protected("results/A04_assembly.bam")
    params:
        cluster_exec=config["cluster_exec"]
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A04_mergeReads.log"
    benchmark: "benchmark/A04_mergeReads.tsv"
    shell:"""
    {{
        echo "##### A04_mergeReads" > {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)
        samtools merge {output.mrg} {input.aln} -@"$numAdditionalThreads"
    }} 2>> {log}
    """

rule A05_indexBam:
    input:
        bam="results/A04_assembly.bam"
    output:
        bai="results/A04_assembly.bam.bai",
        lnk="results/A05_assembly.bai"
    conda: "../envs/sda2.main.yml"
    log: "logs/A05_indexBam.log"
    benchmark: "benchmark/A05_indexBam.tsv"
    params:
        cluster_exec=config["cluster_exec"],
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
    {{
        echo "##### A05_mergeBams" > {log}
        echo "### Index Bam" >> {log}
        samtools index {input.bam} -@{resources.cpus_per_task}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{input.bam}.bai {params.workflowDir}/../{output.lnk}
    }} 2>> {log}
    """
