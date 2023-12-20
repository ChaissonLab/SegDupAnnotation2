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
            bedtools groupby -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14 -c 18,18 -o mean,stdev 1> {output.bed}
        # Out format: chr,start,end,gene,original_chr,original_start,original_end,strand,p_identity,p_accuracy,copy,representative,exons_sizes,exon_starts,depth_by_traditional(non-vcf),depth_by_traditional(non-vcf)-stdDev
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
                {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15/meanCov,$16/meanCov,round($15/meanCov)}}' 1> {output.bed}
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
                {{if ($21==".") \
                    {{$21=0}} \
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$21}}' 1> {output.bed}
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
                ($15>minDepth) \
                    {{print}}' 1> {output.filt}
    }} 2>> {log}
    """

# - F05 Filter by copy number or collapse presence:
#          - autosomes > 2 copies
#          - sex chr > 1 copy
#          - autosomes + sex chr > 1 copy
#          - resolved duplication > 1 copy
#          -> combine with rule MappedSamIdentityDups
rule F05_FindDups:
    input:
        bed="results/F04_resolved_copies_cn2_minDepthFilt.bed",
        sexChrs="results/B04_sex_chrs.txt"
    output:
        dups="results/F05_dups_allFams.bed"
    params:
        sexChrs=expand("{base}",base=config["sex_chr"]), # TODO MAKE SEX CHR SELECTION USER FRIENDLY
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F05_FindDups.log"
    benchmark: "benchmark/F05_FindDups.tsv"
    shell:"""
    {{
        echo "##### F05_FindDups" > {log}
        # sort so OG before copies and all genes adjacent
        # traverse and isolate by gene groups
            # if copyNum of any copy > 1 or 2 keep gene group
            # if resolved copies present keep gene group
        
        # Determine Sex_chrs
        sex_chr_flag=$(cat {input.sexChrs} | tr '\n' ' ' | sed 's/^/-s /g')

        cat {input.bed} | \
            sort -k4,4 -k11,11r -k5,5 -k6,6n -k7,7n -k1,1 -k2,2n -k3,3n | \
            {params.workflowDir}/scripts/F05_SelectDuplications.py /dev/stdin $sex_chr_flag 1> {output.dups}
                     # -s {params.sexChrs} 1> {output.dups}
    }} 2>> {log}
    """

rule F06_GroupIsoformsByNonOverlapping:
    input:
        bed="results/F05_dups_allFams.bed"
    output:
        tsv="results/F06_isoforms_grouped_by_any_overlap.tsv"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F06_GroupIsoformsByNonOverlapping:.log"
    shell:"""
    {{
        echo "##### F06_GroupIsoformsByNonOverlapping:" > {log}
        cat {input.bed} | \
            sort -k1,1 -k2,2n -k3,3n | \
            awk 'BEGIN {{OFS="\t"; chrm=""; groupEnd=0; isos=""; isos_c=""}} \
                (NR==1) \
                    {{chrm=$1; groupEnd=$3}} \
                ($1!=chrm || $2>groupEnd) \
                    {{print isos,isos_c; \
                    chrm=$1; groupEnd=$3; isos=""; isos_c=""}} \
                {{isos=isos $4","; isos_c=isos_c $4"/"$1":"$2"-"$3","; \
                if ($3>groupEnd) {{groupEnd=$3}}}} \
                END \
                    {{print isos,isos_c}}' > {output.tsv} # TODO combine duplicates like last part of NetworkFilter
    }} 2>> {log}
    """

rule F07_PrintRepresentativesOnly:
    input:
        bed="results/F05_dups_allFams.bed"
    output:
        reps="results/F07_dups.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/F07_PrintRepresentativesOnly.log"
    shell:"""
    {{
        echo "##### F07_PrintRepresentativesOnly" > {log}
        cat {input.bed} | \
            awk 'BEGIN {{OFS="\\t"}} ($12=="yes") \
                {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$13,$14,$15,$16,$17,$18}}' 1> {output.reps}
    }} 2>> {log}
    """

