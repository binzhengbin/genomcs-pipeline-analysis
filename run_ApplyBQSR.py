#############################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################

# sample ID
REP_INDEX = {"HS007","HS012","HS013","HS031","HS043","HS066","HS071","HS072","HS080","HS084", \
             "HT011","HT027","HT032","HT044","HT047","HT055","HT056","HT061","HT068","HT077", \
             "HT081","HT089"}
# reference genome pathway
REFERENCE = "/home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/galGal6a/GCA_000002315.5_GRCg6a_genomic.fna "

rule all:
    input:
        expand("6-BQSR/{rep}_bqsr.bam",rep=REP_INDEX),
        expand("6-BQSR/{rep}_bqsr.log",rep=REP_INDEX)

rule ApplyBQSR:
    input:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam",
        "5-baserecalibrator/{rep}_bqsr_data.table"
    output:
        "6-BQSR/{rep}_bqsr.bam"
    params:
        "-Xmx20G"
    log:
        "6-BQSR/{rep}_bqsr.log"
    shell:
        "gatk ApplyBQSR \
        --java-options {params} \
        -R {REFERENCE} \
        -I {input[0]} \
        -O {output} \
        -bqsr {input[1]} \
        > {log} 2>&1 "
