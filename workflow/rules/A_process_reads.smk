# HAPLOTYPE AWARE MODE DISABLED
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

rule A02_faiIndexAsm:
    input:
        asm="results/A01_assembly.fasta"
    output:
        fai="results/A01_assembly.fasta.fai",
        lnk="results/A02_assembly.fai"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_baby,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A02_indexAsm.log"
    benchmark: "benchmark/A02_indexAsm.tsv"
    shell:"""
        echo "##### A02_indexAsm" > {log}
        samtools faidx {input.asm} --fai-idx {output.fai} 2>> {log}

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{output.fai} {params.workflowDir}/../{output.lnk}
    """

rule A03A_mmiIndexAsm:
    input:
        hap="results/A01_assembly.fasta"
    output:
        mmi="results/A03A_assembly.mmi"
    params:
        workflowDir=workflow.basedir,
        read_type=config["read_type"]
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_short"]
    conda: "../envs/sda2.main.yml"
    log: "logs/A03A_mmiIndexAsm.log"
    benchmark: "benchmark/A03A_mmiIndexAsm.tsv"
    shell:"""
        echo "##### A03H_mmiIndexAsm" > {log}
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

rule A03_alignReads:
    input:
        reads=lambda wildcards: bamFiles[wildcards.base],
        mmi=ancient("results/A03A_assembly.mmi")
    output:
        fastq=temp("results/A03_aligned/A03_{base}.fastq"),
        aligned=temp("results/A03_aligned/A03_{base}.bam")
    params:
        read_type=config["read_type"]
    priority: 10
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
        tmp_collate_path=`mktemp -d -p {resources.tmpdir} A03.collate.XXXX.tmp`
        tmp_mm2_path=`mktemp -d -p {resources.tmpdir} A03.mm2.XXXX.tmp`
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} A03.sort.XXXX.tmp`

        if [[ {input.reads} == *".bam" ]]
        then
            echo "### Create fastq files"
            samtools view -h -F 2304 -u -@ "$numAdditionalThreads" {input.reads} | \
                samtools collate -O -u -@ "$numAdditionalThreads" - "$tmp_collate_path" | \
                samtools fastq -@ "$numAdditionalThreads" - > {output.fastq}
        elif [[ {input.reads} == *".fastq.gz" ]]
        then
            echo "### Unzip fastq" >> {log}
            gunzip -c {input.reads} 1> {output.fastq}
        elif [[ {input.reads} == *".fastq" ]]
        then
            cat {input.reads} 1> {output.fastq}
        elif [[ {input.reads} == *".fasta" ]]
        then
            cat {input.reads} 1> {output.fastq}
        elif [[ {input.reads} == *".fa" ]]
        then
            cat {input.reads} 1> {output.fastq}
        fi
        
        echo "### Align Reads" >> {log}
        minimap2 {input.mmi} {output.fastq} -a -x "$mm2_mode" -t {resources.cpus_per_task} --split-prefix "$tmp_mm2_path" | \
            samtools sort -T "$tmp_sort_path"/reads -m "$mem_per_cpu"M -@ "$numAdditionalThreads" -o {output.aligned}

        echo "### Delete Tmp Dirs and their Contents" >> {log}
        rm -rf "$tmp_collate_path" "$tmp_mm2_path" "$tmp_sort_path"
    }} 2>> {log}
    """

rule A04_mergeReads:
    input:
        aln=expand("results/A03_aligned/A03_{base}.bam", base=bamFiles.keys())
    output:
        mrg=protected("results/A04_assembly.bam")
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

rule A05_baiIndexBam:
    input:
        bam="results/A04_assembly.bam"
    output:
        bai="results/A04_assembly.bam.bai",
        lnk="results/A05_assembly.bai"
    conda: "../envs/sda2.main.yml"
    log: "logs/A05_indexBam.log"
    benchmark: "benchmark/A05_indexBam.tsv"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_small,
        cpus_per_task=cluster_cpus_per_task_small,
        runtime=config["cluster_runtime_long"]
    shell:"""
    {{
        echo "##### A05_mergeBams" > {log}
        echo "### Index Bam" >> {log}
        numAdditionalThreads=$(echo "{resources.cpus_per_task} - 1" | bc)
        samtools index {input.bam} -@"$numAdditionalThreads"

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{input.bam}.bai {params.workflowDir}/../{output.lnk}
    }} 2>> {log}
    """
