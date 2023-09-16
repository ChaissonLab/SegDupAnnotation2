# Flow of this smk file:
# - F01 Calculate weighted depth over gene w/o hmm (rule GetGeneCoverage)
# - F02 Label Gene Copy Nums
# - F03 Calculate depth over gene given hmm
# - F04 Filter depth 0.05
# - F05 Filter by copy number or collapse presence: (or at least label these)
#          - autosomal > 2 copies
#          - sex chr > 1 copy
#          -> combine with rule MappedSamIdentityDups
# -     Combine isoforms based on gene name (like combine all VDACs if overlapping like SelectDupsOneIsoform - in addition to E05...)

# Flow for next smk file:
# - Get counts tables
# - Calc summary stats
# - Make figs

rule F01_GetGeneCoverage: # naive depths
    input:
        depths="results/B01_hmm/B01_cov_bins.bed.gz", # presorted
        genes="results/E08_resolved_copies.bed" # presorted
    output:
        bed="results/F01_resolved_copies_cn.bed"
    resources:
        mem_mb=cluster_mem_mb_large,
        cpus_per_task=cluster_cpus_per_task_medium,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/F01_GetGeneCoverage.log"
    benchmark: "benchmark/F01_GetGeneCoverage.tsv"
    shell:"""
    {{
        echo "##### F01_GetGeneCoverage" > {log}
        zcat {input.depths} | \
            bedtools intersect -loj -sorted -a {input.genes} -b /dev/stdin | \
            bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11 -c 15,15 -o mean,stdev 1> {output.bed}
        # Out format: chr,start,end,gene,original_chr,original_start_original_end,strand,p_identity,p_accuracy,copy,depth_by_traditional(non-vcf),depth_by_traditional(non-vcf)-stdDev
    }} 2>> {log}
    """

# Label Copy Nums and convert depth from bp to fraction of mean asm
rule F02_LabelCopyNums:
    input:
        bed="results/F01_resolved_copies_cn.bed",
        mean="results/B03_asm_mean_cov.txt"
    output:
        bed="results/F02_resolved_copies_cn_CNsLabeled.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F02_LabelCopyNums.log"
    benchmark: "benchmark/F02_LabelCopyNums.tsv"
    shell:"""
    {{
        echo "##### F02_LabelCopyNums" > {log}
        cat {input.bed} | \
            awk -v meanCov=$(cat {input.mean}) ' \
                function round(x) {{i=int(x+0.5); \
                                    if (i<1) \
                                        {{return 1}} \
                                    else \
                                        {{return i}}}} \
                BEGIN \
                    {{OFS="\\t"}} \
                {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12/meanCov,$13/meanCov,round($12/meanCov)}}' 1> {output.bed}
    }} 2>> {log}
    """

# TODO verify that "." is 0 and not 1, I think it's 1 but maybe not all the time?
rule F03_GetGeneHmmCoverage: # hmm vcf depths
    input:
        hmm="results/B02_copy_number.bed.gz",
        genes="results/F02_resolved_copies_cn_CNsLabeled.bed" # presorted
    output:
        hmm_noZero=temp("results/F03_copy_number_filtered.bed"),
        bed="results/F03_resolved_copies_cn2.bed"
    resources:
        mem_mb=cluster_mem_mb_large,
        cpus_per_task=cluster_cpus_per_task_medium,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.main.yml"
    log: "logs/F03_GetGeneHmmCoverage.log"
    benchmark: "benchmark/F03_GetGeneHmmCoverage.tsv"
    shell:"""
    {{
        echo "##### F03_GetGeneHmmCoverage" > {log}
        zcat {input.hmm} | \
            awk 'BEGIN {{OFS="\\t"}} ($4!=0) {{print}}' 1> {output.hmm_noZero}
        bedtools intersect -loj -a {input.genes} -b {output.hmm_noZero} -f 1 | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{if ($18==".") \
                    {{$18=0}} \
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$18}}' 1> {output.bed}
        #original rule AddDepthCopyNumber then has bedtools groupby -g 1-4 -c 19 -o max (then cuts and pastes) TODO Delete this line - for reference to original method only
    }} 2>> {log}
    """

rule F04_FilterOutMinDepth:
    input:
        bed="results/F03_resolved_copies_cn2.bed"
    output:
        filt="results/F04_resolved_copies_cn2_minDepthFilt.bed"
    params:
        minDepth=config["min_depth"]
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F04_FilterOutMinDepth.log"
    benchmark: "benchmark/F04_FilterOutMinDepth.tsv"
    shell:"""
    {{
        echo "##### F04_FilterOutMinDepth" > {log}
            cat {input.bed} | \
            awk -v minDepth={params.minDepth} ' \
                BEGIN \
                    {{OFS="\\t"}} \
                ($12>minDepth) \
                    {{print}}' 1> {output.filt}
    }} 2>> {log}
    """

# - F05 Filter by copy number or collapse presence:
#          - autosomes > 2 copies
#          - sex chr > 1 copy
#          - autosomes + sex chr > 1 copy
#          -> combine with rule MappedSamIdentityDups
rule F05_FindDups:
    input:
        bed="results/F04_resolved_copies_cn2_minDepthFilt.bed"
    output:
        dups="results/F05_dups.bed"
    params:
        sexChrs=expand("{base}",base=config["sex_chr"]),
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.python.yml"
    log: "logs/F05_FindDups.log"
    benchmark: "benchmark/F05_FindDups.tsv"
    shell:"""
    {{
        echo "##### F05_FindDups" > {log}
        # sort so OG before copies and all genes adjacent
        # traverse and isolate by gene groups
            # if copyNum of any copy > 1 or 2 keep gene group
            # if resolved copies present keep gene group
        cat {input.bed} | \
            sort -k4,4 -k11,11r -k5,5 -k6,6n -k7,7n -k1,1 -k2,2n -k3,3n | \
            {params.workflowDir}/scripts/F05_SelectDuplications.py /dev/stdin -s {params.sexChrs} 1> {output.dups}
    }} 2>> {log}
    """
