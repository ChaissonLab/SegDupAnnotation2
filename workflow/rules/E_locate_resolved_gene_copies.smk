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
# - E03 filter by identity/accuracy
# - E04 identify exon lengths
# - E05 network approach for overlapping similar sized genes (giving Original priority)
# - E06 paf to bed
# - E07 Annotate Original and sort by gene then by hit loc
# - E08 Filter by Original's hit length margin


rule E01_GetResolvedCopiesPaf:
    input:
        asm="results/A01_assembly.fasta",
        fa="results/D08_resolved_originals.fasta"
    output:
        paf="results/E01_resolved_copies.paf"
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
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_short"],
        tmpdir=tmpDir
    retries: 2
    conda: "../envs/sda2.main.yml"
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

# Locate Exons in gene copies and filter out gene copies that contain <50% of their gene model
rule E04_LocateExons:
    input:
        pafx="results/E03_mapped_resolved_originals_filtered.pafx",
        gm="results/C03_gene_model_filt.fasta",
        asm="results/A01_assembly.fasta"
    output:
        gmidx="results/C03_gene_model_filt.fasta.fai",
        pafxe="results/E04_mapped_resolved_originals_wExons.pafxe"
    params:
        workflowDir=workflow.basedir
    resources:
        mem_mb=cluster_mem_mb_medium,
        cpus_per_task=cluster_cpus_per_task_large,
        runtime=config["cluster_runtime_short"],
        tmpdir=tmpDir
    retries: 2
    conda: "../envs/sda2.main.yml"
    log: "logs/E04_LocateExons.log"
    benchmark: "benchmark/E04_LocateExons.tsv"
    shell:"""
    {{
        echo "##### E04_LocateExons" > {log}
        echo "### Determine Node Variables" >> {log}
        mem_per_cpu="$(echo "{resources.mem_mb}/1.5/{resources.cpus_per_task}" | bc)"
        echo "Memory per cpu: $mem_per_cpu" >> {log}

        echo "### Index Gene Model Fasta" >> {log}
        samtools faidx {input.gm}

        > {resources.tmpdir}/tmp.multiLineTest.err # TODO Delte me

        echo "### Locate Exons" >> {log}
        echo "## Create Locate Function" >> {log}
        locExons () {{ # input = single paf line
            pafLine=`echo "$*" | tr "\\t" " "`

            gene=`echo "$pafLine" | tr "/" "\\t" | awk '(NR==1) {{print $1}}'`
            geneClean=`echo "$pafLine" | tr "/" "\\t" | tr "|" "-" | awk '(NR==1) {{print $1}}'`
            geneLoc=`echo "$pafLine" | tr " " "\\t" | awk '(NR==1) {{print $6":"$8"-"$9}}'`

            tmp_dir=`mktemp -d -p {resources.tmpdir} tmp.getIdent.$geneClean.$$.XXXXXX`

            samtools faidx {input.gm} "$gene" > "$tmp_dir"/gene_model.fasta
            samtools faidx {input.asm} "$geneLoc" > "$tmp_dir"/gene_copy.fasta
            
            minimap2 -x splice -a -t {resources.cpus_per_task} "$tmp_dir"/gene_copy.fasta "$tmp_dir"/gene_model.fasta 1> "$tmp_dir"/gene_copy.sam

            samtools view -b -F 2316 "$tmp_dir"/gene_copy.sam > "$tmp_dir"/gene_copy.bam # -F 4,8,256,2048
            samtools index "$tmp_dir"/gene_copy.bam
            {params.workflowDir}/scripts/D02_FilterMappedLength.py "$tmp_dir"/gene_copy.bam | \
                samtools view -b -o "$tmp_dir"/gene_copy_filt.bam # passes filter if mapped gene copy contains >= 50% gene model

            if [ -f "$tmp_dir"/gene_copy_filt.bam ]; then
                bedtools bamtobed -bed12 -i "$tmp_dir"/gene_copy_filt.bam > "$tmp_dir"/gene_copy.bed
                head -1 "$tmp_dir"/gene_copy.bed | awk -v funcInputPaf="$pafLine" 'BEGIN {{OFS=" "}} {{print funcInputPaf,$11,$12}}' | tr ' ' '\\t' >> {output.pafxe}
                
                cat "$tmp_dir"/gene_copy.bed | tail -n+2 >> {resources.tmpdir}/tmp.multiLineTest.err # TODO Delete me
            fi

            rm -rf "$tmp_dir"
        }}
        export -f locExons

        echo "## Run Locate Function" >> {log}
        cat {input.pafx} | tr "\\t" " " | \
            xargs -P {resources.cpus_per_task} -I % bash -c ' \
                locExons "$@" ' _ %
    }} 2>> {log}
    """

# rule E05_FilterOverlappingGenes: # This rule is now rule D06_FilterOverlappingGenes TODO Delete this block
#     input:
#         pafxe="results/E04_mapped_resolved_originals_wExons.pafxe"
#     output: # TODO should probably sort input first
#         filt="results/E05_mapped_resolved_originals_filtered.pafxe"
#     params:
#         allowOverlappingGenes=config["flag_allow_overlapping_genes"],
#         workflowDir=workflow.basedir
#     localrule: True
#     conda: "../envs/sda2.main.yml"
#     log: "logs/E05_FilterOverlappingGenes.log"
#     benchmark: "benchmark/E05_FilterOverlappingGenes.tsv"
#     shell:"""
#     {{
#         echo "##### E05_FilterOverlappingGenes" > {log}
#         if [ {params.allowOverlappingGenes} = "True" ]
#         then
#             echo "### Do not remove overlapping or intronic genes: create symlink instead." >> {log}
#             ln -s {params.workflowDir}/../{input.pafxe} {params.workflowDir}/../{output.filt}
#         else
#             echo "### Remove intronic and otherwise overlapping genes" >> {log}
#             {params.workflowDir}/scripts/E05_NetworkFilter.py {input.pafxe} | \
#                 awk 'BEGIN {{OFS="\\t"}} ($21==1) {{print $0}}' | \
#                 cut -f1-20 1> {output.filt}
#         fi
#     }} 2>> {log}
#     """

rule E06_FinalResolvedCopiesBed:
    input:
        pafx="results/E04_mapped_resolved_originals_wExons.pafxe"
    output:
        bed="results/E06_mapped_resolved_originals_filtered.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E06_FinalResolvedCopiesBed.log"
    shell:"""
    {{
        echo "##### E06_FinalResolvedCopiesBed" > {log}
        cat {input.pafx} | awk 'BEGIN {{OFS="\\t"}} \
        {{if ($5=="+") \
            {{strand=0}} \
        else \
            {{strand=1}} \
        print $6,$8,$9,$1,strand,$19,$20}}' 1> {output.bed}
        # Out format: chr,start,end,original,strand,p_identity,p_accuracy
    }} 2>> {log}
    """

rule E07_AnnotateOriginal:
    input:
        bed="results/E06_mapped_resolved_originals_filtered.bed"
    output:
        bedAn="results/E07_mapped_resolved_originals_filtered_ogAnnotated.bed"
    params:
        marginToCallOriginal=100
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E07_AnnotateOriginal.log"
    shell:"""
    {{
        echo "##### E07_AnnotateOriginal" > {log}
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

rule E08_FiltByOriginalLengthMargin:
    input:
        bed="results/E07_mapped_resolved_originals_filtered_ogAnnotated.bed"
    output:
        filt="results/E08_resolved_copies.bed"
    params:
        maxLengthMargin=config["max_length_margin"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/E08_FiltByOriginalLengthMargin.log"
    shell:"""
    {{
        echo "##### E08_FiltByOriginalLengthMargin" > {log}
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
