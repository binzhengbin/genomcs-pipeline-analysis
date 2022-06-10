#############################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################

REP_INDEX = {"HS007","HS012","HS013","HS031","HS043","HS066","HS071","HS072","HS080","HS084", \
             "HT011","HT027","HT032","HT044","HT047","HT055","HT056","HT061","HT068","HT077", \
             "HT081","HT089"}
INDEX_BWA = "/home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/galGal6a/galGal6a/"


rule all:
    input:
        expand("fastpresult/afqc_{rep}_1.fq",rep=REP_INDEX),
        expand("fastpresult/afqc_{rep}_2.fq",rep=REP_INDEX),
        expand("fastpresult/unpaired_{rep}_1.fq",rep=REP_INDEX),
        expand("fastpresult/unpaired_{rep}_2.fq",rep=REP_INDEX),
        expand("fastpresult/fastp_{rep}.html",rep=REP_INDEX),
        expand("fastpresult/fastp_{rep}.json",rep=REP_INDEX),
        expand("fastqcresult/{rep}_fastp.log",rep=REP_INDEX)
rule fastp:
    input:
        "RawData/{rep}_1.fq.gz",
        "RawData/{rep}_2.fq.gz"
    output:
        "1-fastpresult/afqc_{rep}_1.fq",
        "1-fastpresult/afqc_{rep}_2.fq",
        "1-fastpresult/unpaired_{rep}_1.fq",
        "1-fastpresult/unpaired_{rep}_2.fq",
        "1-fastpresult/fastp_{rep}.html",
        "1-fastpresult/fastp_{rep}.json"
    log:
        "1-fastqcresult/{rep}_fastp.log"
    shell:
        "fastp --thread 16 --n_base_limit 15 \
        -h {output[4]} -j {output[5]} \
        --qualified_quality_phred 25 \
        --unqualified_percent_limit 50 \
        -i {input[0]} -o {output[0]} \
        -I {input[1]} -O {output[1]} \
        --unpaired1 {output[2]} \
        --unpaired2 {output[3]} \
        > {log} 2>&1 "

rule bwa:
    input:
        "1-fastpresult/afqc_{rep}_1.fq",
        "1-fastpresult/afqc_{rep}_2.fq"
    output:
        "2-bwa_mapping/afqc_{rep}_bwa_gal6a.bam"
    params:
        "@RG\tID:{rep}\tSM:{rep}\tPL:illumina"
    log:
        "2-bwa_mapping/afqc_{rep}_bwa_gal6a.log"
    shell:
        "bwa mem {INDEX_BWA} -t 12 \
        -R {params} \
        {input[0]} {input[1]} | \
        samtools view -S -b - > {output} \
        > {log} 2>&1"

rule bam_file_sort:
    input:
        "2-bwa_mapping/afqc_{rep}_bwa_gal6a.bam"
    output:
        "3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.bam"
    log:
        "3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.log"
    shell:
        "samtolls sort -@ 4 -m 4G -O bam -o {output} {input}"

rule remove_dumplication:
    input:
        "3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.bam"
    output:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam"
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a_metrics.txt"
    log:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.log"
    shell:
        "gatk MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output[0]} \
        --METRICS_FILE {output[1]} \
        --VALIDATION_STRINGENCY SILENT \
        --CREATE_MD5_FILE false \
        --REMOVE_DUPLICATES false \ "

rule Base_Quality_Score_Recalibration:
    input:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam"
    output:
        "5-baserecalibrator/{rep}_bqsr_data.table"
    log:
        "5-baserecalibrator/{rep}_bqsr_data.log"
    shell:
        "gatk BaseRecalibrator \
        -R /home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/galGal6a/GCA_000002315.5_GRCg6a_genomic.fna \
        -I {input} \
        --use-original-qualities \
        -O {output} \
        --known-sites /home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/variation_gallus_gallus/gallus_gallus.20220314.vcf"

rule ApplyBQSR:
    input:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam"
        "5-baserecalibrator/{rep}_bqsr_data.table"
    output:
        "6-BQSR/{rep}_bqsr.bam"
    log:
        "6-BQSR/{rep}_bqsr.log"
    shell:
        "gatk ApplyBQSR \
        -R /home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/galGal6a/GCA_000002315.5_GRCg6a_genomic.fna \
        -I {input[0]} \
        -O {output} \
        -bqsr {input[1]} \
        "
