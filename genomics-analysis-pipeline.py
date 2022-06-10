#############################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################

REP_INDEX = {"007","012","013","031","043","066","071","072","080","084"}
{INDEX_BWA}

rule all:
    input:
        expand("fastpresult/HS_afqc_{rep}_1.fq",rep=REP_INDEX),
        expand("fastpresult/HS_afqc_{rep}_2.fq",rep=REP_INDEX),
        expand("fastpresult/HS_unpaired_{rep}_1.fq",rep=REP_INDEX),
        expand("fastpresult/HS_unpaired_{rep}_2.fq",rep=REP_INDEX),
        expand("fastpresult/HS_fastp_{rep}.html",rep=REP_INDEX),
        expand("fastpresult/HS_fastp_{rep}.json",rep=REP_INDEX),
        expand("fastqcresult/HS_{rep}_fastp.log",rep=REP_INDEX)
rule fastp:
    input:
          "RawData/HS_{rep}_1.fq.gz",
          "RawData/HS_{rep}_2.fq.gz"

    output:
          "fastpresult/HS_afqc_{rep}_1.fq",
          "fastpresult/HS_afqc_{rep}_2.fq",
          "fastpresult/HS_unpaired_{rep}_1.fq",
          "fastpresult/HS_unpaired_{rep}_2.fq",
          "fastpresult/HS_fastp_{rep}.html",
          "fastpresult/HS_fastp_{rep}.json"
    log:
          "fastqcresult/HS_{rep}_fastp.log"
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
        "fastpresult/HS_afqc_{rep}_1.fq",
        "fastpresult/HS_afqc_{rep}_2.fq"
    output:
        "bwa_mapping/HS_afqc_{rep}_bwa_gal6a.bam"
    params:
        "@RG\tID:{rep}\tSM:{rep}\tPL:illumina"
    log:
        "bwa_mapping/HS_afqc_{rep}_bwa_gal6a.log"
    shell:
        "bwa mem {INDEX_BWA} -t 12 \
        -R {params} \
        {input[0]} {input[1]} | \
        samtools view -S -b - > {output} \
        > {log} 2>&1"

rule bam_file_sort:
    input:
        "bwa_mapping/HS_afqc_{rep}_bwa_gal6a.bam"
    output:
        "bwa_sorted/HS_afqc_{rep}_bwa_sorted_gal6a.bam"
    log:
        "bwa_sorted/HS_afqc_{rep}_bwa_sorted_gal6a.log"
    shell:
        "samtolls sort -@ 4 -m 4G -O bam -o {output} {input}"

rule remove_dumplication:
    input:
        "bwa_sorted/HS_afqc_{rep}_bwa_sorted_gal6a.bam"
    output:
        "rmdup/HS_afqc_{rep}_bwa_sorted_rmdup_gal6a.bam"
        "rmdup/HS_afqc_{rep}_bwa_sorted_rmdup_gal6a_metrics.txt"
    params:
        aso = r"queryname"
        so = r"coordinate"
    log:
        "rmdup/HS_afqc_{rep}_bwa_sorted_rmdup_gal6a.log"
    shell:
        "gatk MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output[0]} \
        --METRICS_FILE {output[1]} \
        --VALIDATION_STRINGENCY SILENT \
        --ASSUME_SORT_ORDER '{params.aso}' \
        --CREATE_MD5_FILE false \
        --REMOVE_DUPLICATES false \ "
