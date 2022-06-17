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
        expand("3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.bam",rep=REP_INDEX),
        expand("3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.log",rep=REP_INDEX)


rule bam_file_sort:
    input:
        "2-bwa_mapping/afqc_{rep}_bwa_gal6a.sam"
    output:
        "3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.bam"
    log:
        "3-bwa_sorted/afqc_{rep}_bwa_sorted_gal6a.log"
    shell:
        "samtolls sort -@ 4 -m 20G -O bam -o {output} {input} \
        > {log} 2>&1 "
