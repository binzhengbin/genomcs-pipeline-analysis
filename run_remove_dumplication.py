###########################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################

# sample ID
REP_INDEX = {"HS007","HS012","HS013","HS031","HS043","HS066","HS071","HS072","HS080","HS084", \
             "HT011","HT027","HT032","HT044","HT047","HT055","HT056","HT061","HT068","HT077", \
             "HT081","HT089"}

rule all:
    input:

        expand("4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam",rep=REP_INDEX),
        expand("4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a_metrics.txt",rep=REP_INDEX),
        expand("4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.log",rep=REP_INDEX)

rule remove_dumplication:
    input:
        "3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.bam"
    output:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam",
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a_metrics.txt"
    params:
        "-Xmx20G"
    log:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.log"
    shell:
        "gatk MarkDuplicates \
        --java-options {params} \
        --spark-master local[4] \
        --INPUT {input} \
        --OUTPUT {output[0]} \
        --METRICS_FILE {output[1]} \
        --REMOVE_DUPLICATES false > {log} 2>&1 "
