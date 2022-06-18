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
        expand("7-HaplotypeCaller/{rep}_SNP_only.GCVF.gz",rep=REP_INDEX),
        expand("7-HaplotypeCaller/{rep}_SNP_only.log",rep=REP_INDEX)

rule HaplotypeCaller:
    input:
        "6-BQSR/{rep}_bqsr.bam"
    output:
        "7-HaplotypeCaller/{rep}_SNP_only.GCVF.gz"
    params:
        "-Xmx5G"
    log:
        "7-HaplotypeCaller/{rep}_SNP_only.log"
    shell:
        "gatk HaplotypeCaller \
        --java-options {params} \
        -R {REFERENCE} \
        -ERC GVCF \
        -I {input} \
        -D {SNP_SITES} \
        > {log} 2>&1 "
