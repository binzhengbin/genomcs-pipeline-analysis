#############################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################

# sample ID
REP_INDEX = {"HS007","HS012","HS013","HS031","HS043","HS066","HS071","HS072","HS080","HS084", \
             "HT011","HT027","HT032","HT044","HT047","HT055","HT056","HT061","HT068","HT077", \
             "HT081","HT089"}

# the reference SNP sites from dbsnp
SNP_SITES = "/home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/variation_gallus_gallus/gallus_gallus.20220314.vcf"
# reference genome pathway
REFERENCE = "/home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/galGal6a/GCA_000002315.5_GRCg6a_genomic.fna "

rule all:
    input:
        expand("5-baserecalibrator/{rep}_bqsr_data.table",rep=REP_INDEX),
        expand("5-baserecalibrator/{rep}_bqsr_data.log",rep=REP_INDEX)

rule Base_Quality_Score_Recalibration:
    input:
        "4-rmdup/afqc_{rep}_bwa_sorted_markdup_gal6a.bam"
    output:
        "5-baserecalibrator/{rep}_bqsr_data.table"
    params:
        "-Xmx5G"
    log:
        "5-baserecalibrator/{rep}_bqsr_data.log"
    shell:
        "gatk BaseRecalibratorSpark \
        --java-options {params} \
        --spark-master local[4] \
        -R {REFERENCE} \
        -I {input} \
        -O {output} \
        --known-sites {SNP_SITES} \
        > {log} 2>&1 "
