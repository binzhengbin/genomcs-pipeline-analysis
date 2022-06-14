3301.sh

## build up the pipeline for snp calling

#!/bin/sh

# the reference genome
refGen="/home/RawData/Refseq/Gal6.0/galGal6.fa"  # fasta file required
# the reference SNP sites from dbsnp
snpSites="/home/RawData/Refseq/Gal6.0/gallus_gallus.200123.vcf.gz"

out_dir=$1
sampleN=$2
input1=$3
input2=$4
nt=$5	              # number of cores to be used, specify in the command line

dateR=`date +%F`

# create logfile directory
log_dir=$out_dir"logs/"
mkdir -pm 755 $log_dir

# read group tags
RGgn="YangZhou"       # read group name
RGpl="illumina"       # read group platform
RGsm=$sampleN         # read group sample name

#1. Mapping reads with BWA-MEM, sorting
bam_dir=$out_dir"bamfiles/"
mkdir -pm 755 $bam_dir

timeS1=`date +%s`
bwa mem -aM -t $nt -R "@RG\tID:"$RGgn"\tSM:"$RGsm"\tPL:"$RGpl $refGen $input1 $input2 | \
	samtools sort -m 20G -o $bam_dir$sampleN".sorted.bam" -@ $nt - \
	2>$log_dir$sampleN".bwa."$dateR".log"

timeE1=`date +%s`; echo "step1  Bwa mapping and sorting take $(($timeE1-$timeS1))s"'!'
## when use ";", command1 is independent from command2, both commands will work

#2. Summarizing mapping stats
sum_dir=$out_dir"rawdata_sum/"
mkdir -pm 755 $sum_dir

samtools stats -r $refGen $bam_dir$sampleN".sorted.bam" > $sum_dir$sampleN"sumStats.txt" 
samtools flagstat -@ $nt $bam_dir$sampleN".sorted.bam" > $sum_dir$sampleN"sumFStats.txt" 

timeE2=`date +%s`; echo "step2 Summarizing mapping stats takes $(($timeE2-$timeE1))s"'!'

# 3. Remove Duplicates
RD_dir=$out_dir"remove_duplicates/"
mkdir -pm 755 $RD_dir

/home/software/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx20G" MarkDuplicatesSpark \
     -I $bam_dir$sampleN".sorted.bam" \
     -O $RD_dir$sampleN".RD.bam" \
     -M $RD_dir$sampleN".marked_dup_metrics.txt" \
     2>$log_dir$sampleN".RD_bam."$dateR".log" \
     -- \
     --spark-master local[$nt] 
     #--conf 'spark.executor.cores='$nt \
              
timeE3=`date +%s`; echo "step3 Removing duplicates takes $(($timeE3-$timeE2))s"'!'

# 4a. Base recalibration
RC_dir=$out_dir"Recalibration/"
mkdir -pm 755 $RC_dir

/home/software/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx20G" BaseRecalibratorSpark \
     -I $RD_dir$sampleN".RD.bam" \
	 -R $refGen --known-sites $snpSites \
	 -O $RC_dir$sampleN".RC.table" \
	 2>$log_dir$sampleN".RC_1."$dateR".log" \
	 -- \
     --spark-master local[$nt]
	

# 4b. Base Quality Score Recalibration (BQSR)
/home/software/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx20G" ApplyBQSR \
   	-R $refGen \
   	-I $RD_dir$sampleN".RD.bam" \
   	--bqsr-recal-file $RC_dir$sampleN".RC.table" \
   	-O $RC_dir$sampleN".RC.bam" \
   	2>$log_dir$sampleN".RC_2."$dateR".log" 
    
timeE4=`date +%s` ; echo "step5 BaseRecalibrator takes $(($timeE4-$timeE3))s"'!'

# 5. Variant calling
vcf_dir=$out_dir"vcf/"
mkdir -pm 755 $vcf_dir

/home/software/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx20G" HaplotypeCaller \
	-R $refGen \
	-I $RC_dir$sampleN".RC.bam" \
	-ERC GVCF \
	-D $snpSites \
	-O $vcf_dir$sampleN"_SNP_only.GVCF.gz" \
	2>$log_dir$sampleN".SNPcall_only."$dateR".log"
    
timeE5=`date +%s` ; echo "step6 VariantCalling takes $(($timeE5-$timeE4))s"'!'; echo "The whole process takes $(($timeE5-$timeS1))s"'!'



##02.sh
##Combine the GVCFs
/home/software/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx20G" CombineGVCFs \
   -R /data1/RawData/Refseq/Gal6.0/galGal6.fa \
-V BB1_SNP_only.GVCF.gz \
-V BB2_SNP_only.GVCF.gz \
-V BB3_SNP_only.GVCF.gz \
-O 220502_BLH.GVCF.gz

##Genotype the GVCFs
/data1/software/GATK/gatk4.1.3.sh GenotypeGVCFs --allow-old-rms-mapping-quality-annotation-data \
   -R /data1/RawData/Refseq/Gal6.0/galGal6.fa \
   -V 220605_HT_HS_BLH_TLG.GVCF.gz \
   -O 220608_HT_HS_BLH_TLG.vcf.gz 



##03.sh
#筛选出snp信息
#输入：/home/xjzhang/data_all/merge441.vcf.gz
#输出：/home/xjzhang/data_all/ALL_SNP/merge_SNP.vcf
#/home/xjzhang/data_all/ALL_SNP/01filt_SNP.sh
/home/xjzhang/software/GATK/gatk-4.1.4.1/gatk --java-options "-Xmx40G" SelectVariants \
-V /home/xjzhang/data_all/merge441.vcf.gz \
-select-type SNP \
-O merge_SNP.vcf \
2>merge441_SNP_selectSNP.log


##04.sh
#过滤
less 02merge_SNP.vcf | sed 's/nan/0/g' | bgzip > 02new_merge_SNP.vcf.gz
tabix -p vcf 02new_merge_SNP.vcf.gz

/data1/xjzhang/software/gatk-4.1.4.1/gatk VariantFiltration \
-V 02new_merge_SNP.vcf.gz \
-R /data1/RawData/Refseq/Gal6.0/galGal6.fa \
--filter-expression "QD < 2.0" --filter-name "QD2" \
--filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
--filter-expression "SOR > 3.0" --filter-name "SOR3" \
--filter-expression "FS > 60.0" --filter-name "FS60" \
--filter-expression "MQ < 40.0" --filter-name "MQ40" \
--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O 03new_SNPfilter.vcf.gz \
2>03new_SNPfilter.log

