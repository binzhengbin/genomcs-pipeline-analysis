#############################################################################
# 2022-06-06
# Bin Zheng
# genomics analysis pipeline python script
#############################################################################
# Expecting rule keyword, comment or docstrings inside a rule definition. 这种是因为每行前面有多余的空格，要删除了并一个一个空格输入。


REP_INDEX = {"007","012","013","031","043","066","071","072","080","084"}


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
