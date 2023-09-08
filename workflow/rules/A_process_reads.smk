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
        grid_opts="sbatch -c 1 --mem=32G --time=4:00:00 --partition=qcb --account=mchaisso_100 --output=slurm-logs/slurm-%j.out ", #TODO
        workflowDir=workflow.basedir
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
        grid_opts=config["grid_xlarge"],
        mem_per_cpu_unit="M",
        alt_mem_per_cpu="2000",
        alt_cpus_on_node=1
    resources:
        tmpdir=tmpDir
    conda: "../envs/sda2.main.yml"
    log: "logs/A03_alignReads.{base}.log"
    benchmark: "benchmark/A03_alignReads.{base}.tsv"
    shell:"""
    {{
        echo "##### A03_alignReads" > {log}
        echo "### Determine Node Variables" >> {log}
        if [[ "${{SLURM_CPUS_PER_TASK+defined}}" = defined ]]
        then
            cpus_on_node="$SLURM_CPUS_PER_TASK"
        else
            cpus_on_node={params.alt_cpus_on_node}
        fi
        echo "CPUs on node: $cpus_on_node" >> {log}

        if [[ "${{SLURM_MEM_PER_NODE+defined}}" = defined ]]
        then
            mem_per_cpu="$(echo "($SLURM_MEM_PER_NODE/1.5)/$cpus_on_node" | bc)"
        else
            mem_per_cpu={params.alt_mem_per_cpu}
        fi
        echo "Memory per cpu: $mem_per_cpu" >> {log}

        echo "### Make Tmp Dir" >> {log}
        tmp_sort_path=`mktemp -d -p {resources.tmpdir} reads.sort.XXXX.tmp`

        echo "### Align Reads" >> {log}
        samtools view -h -F 2304 {input.bam} | \
            samtools fastq - | \
            minimap2 {input.asm} - -a -t "$cpus_on_node" | \
            samtools sort -T "$tmp_sort_path"/reads -m "$mem_per_cpu"{params.mem_per_cpu_unit} -@ $(( $cpus_on_node-1 )) -o {output.aligned}

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
        grid_opts=config["grid_medium"],
        alt_cpus_on_node=1
    conda: "../envs/sda2.main.yml"
    log: "logs/A04_mergeReads.log"
    benchmark: "benchmark/A04_mergeReads.tsv"
    shell:"""
    {{
        echo "##### A04_mergeReads" > {log}
        echo "### Determine Node Variables" >> {log}
        if [[ "${{SLURM_CPUS_PER_TASK+defined}}" = defined ]]
        then
            cpus_on_node="$SLURM_CPUS_PER_TASK"
        else
            cpus_on_node={params.alt_cpus_on_node}
        fi
        echo "CPUs on node: $cpus_on_node" >> {log}

        echo "### Merge Bams" >> {log}
        numAdditionalThreads=$(echo "$cpus_on_node - 1" | bc)
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
        grid_opts=config["grid_medium"],
        alt_cpus_on_node=1,
        workflowDir=workflow.basedir
    shell:"""
    {{
        echo "##### A05_mergeBams" > {log}
        echo "### Determine Node Variables" >> {log}
        if [[ "${{SLURM_CPUS_PER_TASK+defined}}" = defined ]]
        then
            cpus_on_node="$SLURM_CPUS_PER_TASK"
        else
            cpus_on_node={params.alt_cpus_on_node}
        fi
        echo "CPUs on node: $cpus_on_node" >> {log}

        echo "### Index Bam" >> {log}
        samtools index {input.bam} -@"$cpus_on_node"

        echo "### Create Link" >> {log}
        ln -s {params.workflowDir}/../{input.bam}.bai {params.workflowDir}/../{output.lnk}
    }} 2>> {log}
    """
