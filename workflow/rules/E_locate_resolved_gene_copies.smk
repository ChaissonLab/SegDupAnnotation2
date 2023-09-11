# Flow of this smk file:
# - D01 map asm to model (getting resolved originals)
# - D02 filter keeping >50% aligned resolved original hits to seq (and index?) (out=bam)
# - D03 bamToBed (out bed)
# - D04 filter out multi exon
# - D05 filter >= MIN_HIT_LENGTH
# - D06 get fasta of resolved original hits
# - D07 properly name resolved original hits fasta
#
#
# - E01 map resolved originals to asm (out: paf)
# - E02 Needleman Wunch alignment of copy hits to original (calc matches/mismatches/etc)
# - E02 calculate resolved copy identity/accuracy
# - E03 filter by identity/accuracy?
# - E04 network approach for overlapping similar sized genes (giving Original priority)
# - E05 paf to bed
# - E06 Annotate Original and sort by gene then by hit loc
# - E07 Filter by Original's hit length margin


rule E01_GetResolvedCopiesPaf:
    input:
        asm="results/A01_assembly.fasta",
        fa="results/D07_resolved_originals.fasta"
    output:
        paf="results/E01_resolved_copies.paf"
    params:
        cluster_exec=config["cluster_exec"]
    resources:
        mem_mb=cluster_mem_mb_xlarge,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/E01_GetResolvedCopiesPaf.log"
    benchmark: "benchmark/E01_GetResolvedCopiesPaf.tsv"
    shell:"""
    {{
        echo "##### E01_GetResolvedCopiesPaf" > {log}
        echo "### Minimap originals to asm" >> {log}
        minimap2 -x asm20 -p 0.2 -N 100 -m 10 -E2,0 -s 10 -t {resources.cpus_per_task} {input.asm} {input.fa} 1> {output.paf}
    }} 2>> {log}
    """

rule E02_GetResolvedCopyIdentities:
    input:
        paf="results/E01_resolved_copies.paf",
        asm="results/A01_assembly.fasta"
    output:
        pafc="results/E02_mapped_resolved_originals.pafxc",
        pafx="results/E02_mapped_resolved_originals.pafx"
    params:
        cluster_exec=config["cluster_exec"],
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_short"],
        tmpdir=tmpDir
    retries: 2
    conda: "../envs/sda2.python.yml"
    log: "logs/E02_GetResolvedCopyIdentities.log"
    benchmark: "benchmark/E02_GetResolvedCopyIdentities.tsv"
    shell:"""
    {{
        echo "##### E02_GetResolvedCopyIdentities" > {log}
        echo "### Align and calculate hit identities" >> {log}
        cat {input.paf} | xargs -P {resources.cpus_per_task} -I % bash -c ' \
            tmp_dir=`mktemp -d -p {resources.tmpdir} tmp.getIdent.$$.XXXXXX`; \
            echo "$@" | tr "\\t" " " > "$tmp_dir"/line.paf; \
            {params.workflowDir}/scripts/E02_CalcPafIdentity.py {input.asm} "$tmp_dir" "$tmp_dir"/line.paf 1>> {output.pafc}; \
            rm -rf "$tmp_dir"; ' _ %
        cat {output.pafc} | cut -f 1-20 1> {output.pafx} # removes cigar string
    }} 2>> {log}
    """

rule E03_FilterLowIdentityPaf:
    input:
        pafx="results/E02_mapped_resolved_originals.pafx"
    output:
        filt="results/E03_mapped_resolved_originals_filtered.pafx"
    params:
        min_copy_identity=config["min_copy_identity"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E03_FilterLowIdentityPaf.log"
    benchmark: "benchmark/E03_FilterLowIdentityPaf.tsv"
    shell:"""
    {{
        echo "##### E03_FilterLowIdentityPaf" > {log}
        cat {input.pafx} | \
            awk -v minId={params.min_copy_identity} ' \
                BEGIN \
                    {{OFS="\\t"}} \
                ($19 >= minId) \
                    {{print $0}}' | \
            sort -k1,1 -k6,6 -k8,8n -k9,9n 1> {output.filt}
    }} 2>> {log}
    """

rule E04_FilterOverlappingGenes:
    input:
        pafx="results/E03_mapped_resolved_originals_filtered.pafx"
    output: # TODO should probably sort input first
        filt="results/E04_mapped_resolved_originals_filtered.pafx"
    params:
        allowOverlappingGenes=config["flag_allow_overlapping_genes"],
        workflowDir=workflow.basedir
    conda: "../envs/sda2.python.yml"
    localrule: True # TODO Double check
    log: "logs/E04_FilterOverlappingGenes.log"
    benchmark: "benchmark/E04_FilterOverlappingGenes.tsv"
    shell:"""
    {{
        echo "##### E04_FilterOverlappingGenes" > {log}
        if [ {params.allowOverlappingGenes} = "True" ]
        then
            echo "### Do not remove overlapping or intronic genes: create symlink instead." >> {log}
            ln -s {params.workflowDir}/../{input.pafx} {params.workflowDir}/../{output.filt}
        else
            echo "### Remove intronic and otherwise overlapping genes" >> {log}
            cat {input.pafx} | \
                awk 'BEGIN {{OFS="\\t"}} {{print $0,NR,"0"}}' | \
                {params.workflowDir}/scripts/E04_NetworkFilter.py /dev/stdin | \
                awk 'BEGIN {{OFS="\\t"}} ($22==1) {{print $0}}' | \
                cut -f1-20 1> {output.filt} # TODO Pick and tune optimal networking alg
        fi
    }} 2>> {log}
    """

rule E05_FinalResolvedCopiesBed:
    input:
        pafx="results/E04_mapped_resolved_originals_filtered.pafx"
    output:
        bed="results/E05_mapped_resolved_originals_filtered.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E05_FinalResolvedCopiesBed.log"
    shell:"""
    {{
        echo "##### E05_FinalResolvedCopiesBed" > {log}
        cat {input.pafx} | awk 'BEGIN {{OFS="\\t"}} \
        {{if ($5=="+") \
            {{strand=0}} \
        else \
            {{strand=1}} \
        print $6,$8,$9,$1,strand,$19,$20}}' 1> {output.bed}
        # Out format: chr,start,end,original,strand,p_identity,p_accuracy
    }} 2>> {log}
    """

rule E06_AnnotateOriginal:
    input:
        bed="results/E05_mapped_resolved_originals_filtered.bed"
    output:
        bedAn="results/E06_mapped_resolved_originals_filtered_ogAnnotated.bed"
    params:
        marginToCallOriginal=100
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E06_AnnotateOriginal.log"
    shell:"""
    {{
        echo "##### E06_AnnotateOriginal" > {log}
        cat {input.bed} | \
            tr '/' '\\t' | \
            awk -v margin={params.marginToCallOriginal} ' \
                function abs(v) {{v += 0; return v < 0 ? -v : v}} \
                BEGIN \
                    {{OFS="\\t"}} \
                {{split($5,hit,":"); \
                split(hit[2],hitLocs,"-"); \
                og_chr=$1; og_start=$2; og_end=$3; \
                hit_chr=hit[1]; hit_start=hitLocs[1]; hit_end=hitLocs[2]; \
                if (og_chr==hit_chr && abs(og_start-hit_start)<margin && abs(og_end-hit_end)<margin) \
                    {{type="Original"}} \
                else \
                    {{type="Copy"}} \
                print $1,$2,$3,$4,hit_chr,hit_start,hit_end,$6,$7,$8,type}}' | \
                sort -k1,1 -k2,2n -k3,3n -k4,4 1> {output.bedAn}
        # Out format: chr,start,end,gene,original_chr,original_start_original_end,strand,p_identity,p_accuracy,copy # sorted by chrom, start, end, then gene
    }} 2>> {log}
    """

rule E07_FiltByOriginalLengthMargin:
    input:
        bed="results/E06_mapped_resolved_originals_filtered_ogAnnotated.bed"
    output:
        filt="results/E07_resolved_copies.bed"
    params:
        maxLengthMargin=config["max_length_margin"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E07_FiltByOriginalLengthMargin.log"
    shell:"""
    {{
        echo "##### E07_FiltByOriginalLengthMargin" > {log}
        cat {input.bed} | \
            awk -v maxLengthMargin={params.maxLengthMargin} ' \
                BEGIN \
                    {{OFS="\\t"; \
                    highMarg=1+maxLengthMargin; \
                    lowMarg=1-maxLengthMargin}} \
                (($3-$2)<highMarg*($7-$6) && ($3-$2)>lowMarg*($7-$6)) \
                    {{print $0}}' | \
                sort -k1,1 -k2,2n -k3,3n -k4,4 1> {output.filt}
        # Out format: chr,start,end,gene,original_chr,original_start_original_end,strand,p_identity,p_accuracy,copy # sorted by chrom, start, end, then gene
    }} 2>> {log}
    """
