# Flow of this smk file:
# - G01 Add header to dup results
# -     Get counts tables
# -     Calc summary stats
# -     Make figs

rule G01_AddHeader:
    input:
        dups="results/F07_dups.bed"
    output:
        bed="results/G01_dups.bed",
        bed12="results/G01_dups.igv.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/G01_AddHeader.log"
    shell:"""
    {{
        echo "##### G01_AddHeader" > {log}
        cat {input.dups} | \
            sort -k4,4 -k11,11r -k5,5 -k6,6n -k7,7n -k1,1 -k2,2n -k3,3n | \
            awk ' \
                BEGIN \
                    {{OFS="\\t"; \
                    print "#chr","start","end","gene","orig_chr","orig_start","orig_end", \
                    "strand","p_identity","p_accuracy","identity", \
                    "depth","depth_stdev","copy_num","depth_by_vcf"}} \
                {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$14,$15,$16,$17}}' 1> {output.bed}

        echo "### IGV BED File" >> {log}
        cat {input.dups} | \
            awk 'BEGIN {{OFS="\\t"}} \
                {{split($12,exonSizes,","); \
                print $1,$2,$3,$4,int($9*1000),$8,$2,$3,"255,0,0",length(exonSizes),$12,$13}}' 1> {output.bed12}
    }} 2>> {log}
    """

# counting resolved copies of genes and collapses of genes EXCLUDING ORIGINALS
# Note one copy per line, but Originals (or if Original filtered out,
# first copy removed - so it only notes high quality copies.)
rule G02_GeneCountFact:
    input:
        dups="results/G01_dups.bed"
    output:
        fact="results/G02_dups_fact.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/G02_GeneCountFact.log"
    shell:"""
    {{
        echo "##### G02_GeneCountFact" > {log}
        cat {input.dups} | \
            awk 'BEGIN \
                    {{OFS="\\t"; \
                    lastGene=""}} \
                (NR==1) \
                    {{print $0,"resolved"}} \
                (NR>1) \
                    {{if (lastGene==$4) \
                        {{print $0,"multi"}} \
                    nCopy=int($14); \
                    for (i=1; i<nCopy; i++) \
                        {{print $0,"collapse"}} \
                    lastGene=$4}}' 1> {output.fact}
    }} 2>> {log}
    """

rule G03_PerGeneCounts:
    input:
        dups="results/G01_dups.bed",
        mean="results/B03_asm_mean_cov.txt"
    output:
        genes="results/G03_per_gene_counts.bed"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/G03_PerGeneCounts.log"
    shell:"""
    {{
        echo "##### G03_PerGeneCounts" > {log}
        # Assumption: input data sorted by gene first
        cat {input.dups} | \
            awk -v meanCov=$(cat {input.mean}) ' \
                BEGIN \
                    {{OFS="\\t"; \
                    gene=""; depth_sum=0; count=0; resolvedNum=0;}} \
                (NR==1) \
                    {{print "gene","depthNoramlizedToAsm","depth","measuredCov","copyCount","resolvedCount"}} \
                (NR==2) \
                    {{gene=$4; depth_sum=$12; count=$14; resolvedNum=1}} \
                (NR>2) \
                    {{if ($4==gene) \
                        {{depth_sum+=$12; count+=$14; resolvedNum++}} \
                    else \
                        {{print gene,depth_sum,depth_sum*meanCov,sprintf("%1.f",depth_sum),count,resolvedNum; \
                        gene=$4; depth_sum=$12; count=$14; resolvedNum=1}} }} \
                END \
                    {{print gene,depth_sum,depth_sum*meanCov,sprintf("%1.f",depth_sum),count,resolvedNum;}}' 1> {output.genes}
    }} 2>> {log}
    """

rule G04_SummaryStats:
    input:
        fact="results/G02_dups_fact.bed",
        mean="results/B03_asm_mean_cov.txt",
        vcfDups="results/B02_copy_number.bed.gz"
    output:
        summary="results/G04_summary_stats.tsv"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/G04_SummaryStats.log"
    benchmark: "benchmark/G04_SummaryStats.tsv"
    shell:"""
    {{
        echo "##### G04_SummaryStats" > {log}
        collapsedBasesTotal="$(zcat {input.vcfDups} | awk \
            'BEGIN {{OFS="\\t";sum=0}} \
            ($4>1) \
                {{sum+=($3-$2)*($4-1)}} \
            END {{print sum}}')" # sum of total collapsed bases (regardless of occurence within gene - calculated using hmm vcf file)
        collapsedBasesInGenes="$(cat {input.fact} | awk \
            'BEGIN {{OFS="\\t"; sum=0}} \
            (NR>1 && $16=="collapse") \
                {{sum+=($3-$2)*$12}} \
            END {{printf "%1.0i\\n", sum}}')" # sum of collapsed gene copies (excludes 'original' copies and resolved copies)
        resolvedBasesinCollapsedGenes="$(cat {input.fact} | awk \
            'BEGIN {{OFS="\\t"; sum=0}} \
            (NR>1 && $16=="collapse") \
                {{sum+=($3-$2)}} \
            END {{printf "%1.0i\\n", sum}}')" # sum of collapsed gene copies (excludes 'original' copies and resolved copies)
        resolvedBases="$(cat {input.fact} | awk \
            'BEGIN {{OFS="\\t"; sum=0}} \
            (NR>1 && $16=="multi") \
                {{sum+=($3-$2)}} \
            END {{print sum}}')" # sum of resolved copies of genes (excludes 'original' copies and collapses)
        resolvedGeneCount="$(cat {input.fact} | awk \
            'BEGIN {{OFS="\\t"}} \
            (NR>1 && $16=="multi") \
                {{print $4}}' | \
            uniq | wc -l)"
        collapsedGeneCount="$(cat {input.fact} | awk \
            'BEGIN {{OFS="\\t"}} \
            (NR>1 && $16=="collapse") \
                {{print $4}}' | \
            uniq | wc -l)"
        collapsedDuplicationCount="$(cat {input.fact} | awk \
            'BEGIN {{OFS="\\t"; count=0}} \
            (NR>1 && $16=="collapse") \
                {{count+=1}} \
            END {{print count}}')"

        echo -e "resolved_duplicated_gene_count\\t$resolvedGeneCount\\n\
resolved_duplicate_bases\\t$resolvedBases\\n\
collapsed_duplications_count\\t$collapsedDuplicationCount\\n\
collapsed_duplicated_gene_count\\t$collapsedGeneCount\\n\
bases_in_collapsed_genes\\t$collapsedBasesInGenes\\n\
resolved_bases_in_collapsed_genes\\t$resolvedBasesinCollapsedGenes\\n\
Total_collapsed_bases\\t$collapsedBasesTotal" > {output.summary}
    }} 2>> {log}
    """

rule G05_SummaryFigs:
    input:
        dups="results/G01_dups.bed",
        mean="results/B03_asm_mean_cov.txt",
    output:
        comboTmp=temp("results/G05_dups.tsv"),
        resTmp=temp("results/G05_dups_resolved.tsv"),
        colTmp=temp("results/G05_dups_collapsed.tsv"),
        plot_combo="results/G05_depthPlot_combo.pdf",
        plot_res="results/G05_depthPlot_res.pdf",
        plot_col="results/G05_depthPlot_col.pdf",
        plot_merged="results/G05_depthPlot_merged.pdf",
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.r.yml"
    log: "logs/G05_SummaryFigs.log"
    benchmark: "benchmark/G05_SummaryFigs.tsv"
    shell:"""
    {{
        echo "##### G05_SummaryFigs" > {log}
        echo "### Reformat inputs" >> {log}
        # Below assumes sorted by gene name
        cat {input.dups} | \
            awk 'BEGIN {{OFS="\\t"; gene=""; gene_num=0}} \
                (NR==1) \
                    {{print "gene_num",$0}} \
                (NR>1) \
                    {{if (gene!=$4) \
                        {{geneNum+=1}} \
                    gene=$4; \
                    print geneNum,$0}}' | \
            tr -d '#' > {output.comboTmp}
        cat {input.dups} | \
            awk 'BEGIN {{OFS="\\t"; gene=""; gene_num=0}} \
                (NR==1) \
                    {{print "gene_num",$0}} \
                (NR>1 && $13<=1) \
                    {{if (gene!=$4) \
                        {{geneNum+=1}} \
                    gene=$4; \
                    print geneNum,$0}}' | \
            tr -d '#' > {output.resTmp}
        cat {input.dups} | \
            awk 'BEGIN {{OFS="\\t"; gene=""; gene_num=0}} \
                (NR==1) \
                    {{print "gene_num",$0}} \
                (NR>1 && $13>1) \
                    {{if (gene!=$4) \
                        {{geneNum+=1}} \
                    gene=$4; \
                    print geneNum,$0}}' | \
            tr -d '#' > {output.colTmp}

        echo "### Generate figs" >> {log}
        Rscript {params.workflowDir}/scripts/G05_SummaryFigs.R {params.workflowDir} {output.comboTmp} {output.resTmp} {output.colTmp}
    }} 2>> {log}
    """