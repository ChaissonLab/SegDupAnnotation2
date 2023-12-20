rule B01_RunHmm:
    input:
        asm="results/A01_assembly.fasta",
        fai="results/A01_assembly.fasta.fai",
        bam=ancient("results/A04_assembly.bam"),
        bai="results/A05_assembly.bai" # TODO necessary?
    output:
       vcf="results/B01_hmm/B01_copy_number.vcf",
       cov=temp("results/B01_hmm/B01_cov_bins.bed"),
       covz="results/B01_hmm/B01_cov_bins.bed.gz",
       tbi="results/B01_hmm/B01_cov_bins.bed.gz.tbi", # TODO is this useful?
       snvs="results/B01_hmm/B01_snvs.tsv"
    resources:
        mem_mb=cluster_mem_mb_large,
        cpus_per_task=cluster_cpus_per_task_medium,
        runtime=config["cluster_runtime_long"]
    conda: "../envs/sda2.hmcnc.yml"
    log: "logs/B01_RunHmm.log"
    benchmark: "benchmark/B01_RunHmm.tsv"
    shell:"""
    {{
        echo "##### B01_RunHmm" > {log}
        echo "### Run HMM" >> {log}
        hmmcnc {input.asm} -a {input.bam} -t {resources.cpus_per_task} -B {output.cov} -S {output.snvs} -o {output.vcf}

        echo "### Sort and Index Output" >> {log}
        cat {output.cov} | \
            sort -k1,1 -k2,2n -k3,3n | \
            bgzip -c 1> {output.covz}
        tabix {output.covz}
    }} 2>> {log}
    """

rule B02_GetCN:
    input:
       vcf="results/B01_hmm/B01_copy_number.vcf"
    output:
       bed="results/B02_copy_number.bed.gz"
    params:
        workflowDir=workflow.basedir
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/B02_GetCN.log"
    benchmark: "benchmark/B02_GetCN.tsv"
    shell:"""
    {{
        echo "##### B02_GetCN" > {log}
        {params.workflowDir}/scripts/B02_CovVcfToBed.py {input.vcf} | \
            gzip -c 1> {output.bed}
        # Out format: chr,start,end,copy_num,read_depth
    }} 2>> {log}
    """

rule B03_CalcMeanCov:
    input:
        cov="results/B01_hmm/B01_cov_bins.bed.gz"
    output:
        mc="results/B03_asm_mean_cov.txt"
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/B03_CalcMeanCov.log"
    benchmark: "benchmark/B03_CalcMeanCov.tsv"
    shell:"""
    {{
        echo "##### B03_CalcMeanCov" > {log}
        zcat {input.cov} | \
            awk '{{sum+=$4; count+=1}} \
                END {{print sum/count}}' 1> {output.mc}
    }} 2>> {log}
    """

rule B04_CalcMeanDepthPerChrom:
    input:
        cov="results/B01_hmm/B01_cov_bins.bed.gz",
        mean="results/B03_asm_mean_cov.txt"
    output:
        tsv="results/B04_chr_mean_cov.tsv",
        sexChrs="results/B04_sex_chrs.txt"
    params:
        min_sex_chr_len=10000000, # =10Mb
        depth_margin=8 # TODO CHANGE THIS, Make it proportional or something
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/B04_CalcMeanChrPerChrom.log"
    benchmark: "benchmark/B04_CalcMeanChrPerChrom.tsv"
    shell:"""
    {{
        echo "##### B04_CalcMeanChrPerChrom" > {log}
        # Assumes bin file sorted by chr name
        zcat {input.cov} | \
            awk 'BEGIN \
                    {{OFS="\t"; \
                    chr=""; depth_sum=0; count=0; chrDepth=0; \
                    print "chr","mean_depth","length"}} \
                (NR==1) \
                    {{chr=$1}} \
                (chr!=$1) \
                    {{chrDepth=depth_sum/count; \
                    print chr,chrDepth,count*100; \
                    chr=$1; depth_sum=0; count=0}} \
                (chr==$1) \
                    {{depth_sum+=$4; count++}} \
                END \
                    {{print chr,chrDepth,count*100}}' > {output.tsv}

        cat {output.tsv} | \
            tail -n+2 | \
            awk -v minLen={params.min_sex_chr_len} -v meanCov=$(cat {input.mean}) -v depthMargin={params.depth_margin} \
                'BEGIN {{OFS="\\t"}} \
                ($3>=minLen && $2<meanCov/2+depthMargin && $2>meanCov/2-depthMargin) \
                    {{print $1}}' 1> {output.sexChrs}
    }} 2>> {log}
    """

# TODO Delete below, DO NOT USE
# rule B03_CombineCNRanges:
#     input:
#        bed="results/B02_copy_number.bed.gz"
#     output:
#        comb="results/B03_hmm_copy_number.bed.gz",
#        test=touch("B.done")
#     localrule: True
#     log: "logs/B03_CombineCNRanges.log"
#     benchmark: "benchmark/B03_CombineCNRanges.tsv"
#     shell:"""
#         echo "##### B03_CombineCNRanges" > {log}
#         zcat {input.bed} | \
#         bedtools merge -i /dev/stdin -c 4,4,4 -o min,max,mean | \
#             gzip -c 1> {output.comb} 2>> {log} # really bad method should be a weighted mean!!
#    """
