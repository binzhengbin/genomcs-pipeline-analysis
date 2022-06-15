#############################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################

# sample ID
REP_INDEX = {"HS007","HS012","HS013","HS031","HS043","HS066","HS071","HS072","HS080","HS084", \
             "HT011","HT027","HT032","HT044","HT047","HT055","HT056","HT061","HT068","HT077", \
             "HT081","HT089"}
# bwa index pathway
INDEX_BWA = "/home/bzheng/genomics-data/X101SC21110256-Z01-J031/ref/galGal6a/galGal6a"

rule all:
    input:
        expand("2-bwa_mapping/afqc_{rep}_bwa_gal6a.sam",rep=REP_INDEX),
        expand("2-bwa_mapping/afqc_{rep}_bwa_gal6a.log",rep=REP_INDEX),

rule bwa:
    input:
        "1-fastpresult/afqc_{rep}_1.fq",
        "1-fastpresult/afqc_{rep}_2.fq"
    output:
        "2-bwa_mapping/afqc_{rep}_bwa_gal6a.sam"
    params:
        "@RG\tID:XinHua\tSM:{rep}\tPL:illumina"
    log:
        "2-bwa_mapping/afqc_{rep}_bwa_gal6a.log"
    shell:
        "bwa mem {INDEX_BWA} -t 12 \
        -aM \
        -R {params} \
        {input[0]} {input[1]} > {output} \
        > {log} 2>&1"
