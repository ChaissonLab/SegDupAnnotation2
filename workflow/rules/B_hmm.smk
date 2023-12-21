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
        cov="results/B01_hmm/B01_cov_bins.bed.gz"
    output:
        tsv="results/B04_chr_mean_cov.tsv"
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
                    print "chr","depth_in_bp","length"}} \
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
    }} 2>> {log}
    """

rule B05_IdentifyHaploidChrs:
    input:
        tsv="results/B04_chr_mean_cov.tsv",
        mean="results/B03_asm_mean_cov.txt"
    output:
        hapChrs="results/B05_haploid_chrs.txt"
    params:
        autodetect_haploid_chrs_flag=config["flag_autodetect_haploid_chrs"],
        min_sex_chr_len=10000000, # =10Mb
        depth_margin=10, # as percentage of normalized depth
        hapChrs=expand("{base}",base=config["haploid_chrs"])
    localrule: True
    conda: "../envs/sda2.main.yml"
    log: "logs/B05_IdentifyHaploidChrs.log"
    shell:"""
    {{
        echo "##### B05_IdentifyHaploidChrs" > {log}
        if [ {params.autodetect_haploid_chrs_flag} = "True" ]
        then
            echo "### Autodetecting haploid chromosomes" >> {log}
            cat {input.tsv} | \
                tail -n+2 | \
                awk -v minLen={params.min_sex_chr_len} -v meanCov=$(cat {input.mean}) -v depthMargin={params.depth_margin} \
                    'BEGIN {{OFS="\\t"}} \
                    ($3>=minLen && $2/meanCov*100<50+depthMargin && $2/meanCov*100>50-depthMargin) \
                        {{print $1}}' 1> {output.hapChrs}
        else
            echo "### Using supplied haploid chromosomes" >> {log}
            echo "# If no haploid chromosomes were supplied, all chromosomes assumed to be diploid."
            echo {params.hapChrs} | tr ' ' '\n' 1> {output.hapChrs}
        fi
    }} 2>> {log}
    """
    