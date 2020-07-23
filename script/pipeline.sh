# nohup ./sample_fq_to_bam.sh /work_dir/ sample-id > log.txt 2 > &1 &
# command-line arguments, my work_dir = /home/ruizhi/data/wgs_pipeline/
# fastq raw data is stored at ${work_dir}/fasta_file
# assume your fastq filename format is: SAMPLE-ID_*.1.fastq.gz   SAMPLE-ID_*.2.fastq.gz
work_dir=$1
sample=$2

#reference
reference=${work_dir}/reference/genome/hg38.fa
GATK_bundle=${work_dir}/resources/gatk/bundle/hg38
# create result folders
if [ ! -d ${work_dir}/result/qc/${sample}_result/ ]
then mkdir -p ${work_dir}/result/qc/${sample}_result
fi

if [ ! -d ${work_dir}/result/qc/${sample}_result/ ]
then mkdir -p ${work_dir}/result/qc/${sample}_result
fi

if [ ! -d ${work_dir}/result/cleanfq/${sample}_result/ ]
then mkdir -p ${work_dir}/result/cleanfq/${sample}_result/
fi

if [ ! -d ${work_dir}/result/bwa/${sample}_result/ ]
then mkdir -p ${work_dir}/result/bwa/${sample}_result/
fi

if [ ! -d ${work_dir}/result/bwa/${sample}_result/tmp ]
then mkdir -p ${work_dir}/result/bwa/${sample}_result/tmp
fi

if [ ! -d ${work_dir}/result/gatk ]
then mkdir -p ${work_dir}/result/gatk
fi

if [ ! -d ${work_dir}/result/gatk/${sample}_result/VCF ]
then mkdir -p ${work_dir}/result/gatk/${sample}_result/VCF
fi

if [ ! -d ${work_dir}/result/gatk/${sample}_result/tmp ]
then mkdir -p ${work_dir}/result/gatk/${sample}_result/tmp
fi

if [ ! -d ${work_dir}/result/gatk ]
then mkdir -p ${work_dir}/result/gatk
fi

if [ ! -d ${work_dir}/result/gatk/${sample}_result/VCF ]
then mkdir -p ${work_dir}/result/gatk/${sample}_result/VCF
fi

if [ ! -d ${work_dir}/result/gatk/${sample}_result/tmp ]
then mkdir -p ${work_dir}/result/gatk/${sample}_result/tmp
fi
echo "umi_tools extract sarted at $(date)" 
cd /mnt/lab-in/whole-genome-sequencing/103734-002/
for read1 in *${sample}*R1.fastq.gz;do
    read2=$(echo $read1| sed 's/R1.fastq.gz/R2.fastq.gz/');
    read3=$(echo $read1| sed 's/R1.fastq.gz/R3.fastq.gz/');
    output1=$(echo $read1| sed 's/R1.fastq.gz/R1_umi/');
    output2=$(echo $read3| sed 's/R3.fastq.gz/R2_umi/');
    umi_tools extract --bc-pattern=NNNNNNN --stdin=$read2 --read2-in=$read1 --stdout=${work_dir}/fasta_file/${output1}.fastq.gz --read2-stdout;
    umi_tools extract --bc-pattern=NNNNNNN --stdin=$read2 --read2-in=$read3 --stdout=${work_dir}/fasta_file/${output2}.fastq.gz --read2-stdout;   
done

cat ${work_dir}/fasta_file/*${sample}*R1_umi.fastq.gz > ${work_dir}/fasta_file/${sample}_R1.fastq.gz
cat ${work_dir}/fasta_file/*${sample}*R2_umi.fastq.gz > ${work_dir}/fasta_file/${sample}_R2.fastq.gz
echo "umi_tools extract ended at $(date)" 

rm ${work_dir}/fasta_file/*${sample}*R1_umi.fastq.gz ${work_dir}/fasta_file/*${sample}*R2_umi.fastq.gz 


# fastqc: check fastq files quality, result stored at result/qc/${sample}_result/
echo "fastqc checking data quality started at $(date)" 
cd ${work_dir}/result/qc/${sample}_result/
fastqc -t 30 ${work_dir}/fasta_file/${sample}*gz -o ${work_dir}/result/qc/${sample}_result/
echo "fastqc checking data quality finished at $(date)" 

# filter low_quality data and remove adapter
cd ${work_dir}/result/cleanfq/${sample}_result/
echo "trim_galore cutting adaptors started at $(date)" 
ls ${work_dir}/fasta_file/${sample}*1.*.gz > ./1;cat 1
ls ${work_dir}/fasta_file/${sample}*2.*.gz > ./2;cat 2
paste 1 2 > config1 && rm 1 2
cat ${work_dir}/result/cleanfq/${sample}_result/config1 | while read id; do
   arr=(${id}); fq1=${arr[0]}; fq2=${arr[1]}; trim_galore -j 20 --fastqc -q 20 --phred33 --length 36 -e 0.1 --stringency 3 --paired $fq1 $fq2 -o ${work_dir}/result/cleanfq/${sample}_result/
done
echo "trim_galore cutting adapters finished at $(date)" && rm config1

# bwa mem: mapping reads to GRCh38
echo "bwa mem mapping reads started at $(date)"
cd ${work_dir}/result/bwa/${sample}_result/
ls ${work_dir}/result/cleanfq/${sample}_result/${sample}_R1_val_1.fq.gz > ./1;cat 1
ls ${work_dir}/result/cleanfq/${sample}_result/${sample}_R2_val_2.fq.gz > ./2;cat 2
paste 1 2 > config2 && rm 1 2
cat ${work_dir}/result/bwa/${sample}_result/config2 | while read id; do
    arr=(${id}); fq1=${arr[0]}; fq2=${arr[1]}; tmp=`basename ${arr[0]}`; sample=${sample}; bwa mem -M -v 2 -t 20 -R "@RG\tID:${sample}\tSM:${sample}\tLB:WGS\tPL:Illumina" $reference $fq1 $fq2 | samtools view -@ 30 -bS - > ${sample}.bam
done
echo "bwa mapping reads finished at $(date)" && rm config2

echo "sambamba flagstat started at $(date)"
sambamba flagstat -t 20 ${sample}.bam > ${sample}.flagstat.tsv
echo "sambamba flagstat finished at $(date)" 

# sort bam files
echo "sambamba sort sorting bam files started at $(date)"
sambamba sort -t 20 -m 2G --tmpdir=${work_dir}/result/bwa/${sample}_result/tmp -o ${work_dir}/result/bwa/${sample}_result/${sample}.sorted.bam ${sample}.bam
echo "sambamba sort sorting bam files finished at $(date)"

#mark duplicated sequences, it'll generate bai index files automatically
echo "umi_tools deduplication marks duplicated sequences started at $(date)"
cd ${work_dir}/result/bwa/${sample}_result/
umi_tools dedup -I ${sample}.sorted.bam -S ${work_dir}/result/gatk/${sample}_result/${sample}.sorted.markdup.bam 
echo "umi_tools deduplication marks duplicated sequences started at $(date)"
 
# excute BQSR
cd ${work_dir}/result/gatk/${sample}_result/
echo "gatk BaseRecalibrator recalibrating started at $(date)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" BaseRecalibrator \
    -R $GATK_bundle/Homo_sapiens_assembly38.fasta \
    -I ${sample}.sorted.markdup.bam  \
    --tmp-dir tmp \
    --known-sites $GATK_bundle/dbsnp_146.hg38.vcf.gz \
    --known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O ${sample}.sorted.markdup.recal_data.table
echo "gatk BaseRecalibrator recalibrating finished at $(date)"

# gatk ApplyBQSR 
echo "gatk ApplyBQSR started at $(date)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" ApplyBQSR \
   --bqsr-recal-file ${sample}.sorted.markdup.recal_data.table \
   --tmp-dir tmp \
   -R $GATK_bundle/Homo_sapiens_assembly38.fasta \
   -I ${sample}.sorted.markdup.bam \
   -O ${sample}.sorted.markdup.BQSR.bam
echo "gatk ApplyBQSR finished at $(date)"

cd ${work_dir}/result/gatk/${sample}_result/
# generate raw VCF files
echo "gatk HaplotypeCaller generating VCF started at $(date)"
chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
for a in ${sample}; do
    for i in ${chrom[@]}; do
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" HaplotypeCaller \
            --emit-ref-confidence GVCF \
            -R $GATK_bundle/Homo_sapiens_assembly38.fasta  \
            --tmp-dir tmp \
            -I ${sample}.sorted.markdup.BQSR.bam \
            --native-pair-hmm-threads 8 \
            -L $i \
            -O ${sample}.HC.${i}.g.vcf.gz && \
        gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" GenotypeGVCFs \
            -R $GATK_bundle/Homo_sapiens_assembly38.fasta  \
            -V ${sample}.HC.${i}.g.vcf.gz \
            -O ${sample}.HC.${i}.vcf.gz &
    done
    wait
done

for i in ${chrom[@]}; do
    ls ${sample}.HC.${i}.vcf.gz >> ${sample}.input_variant_files.list
done

gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" MergeVcfs \
    --TMP_DIR=tmp \
    -I=${sample}.input_variant_files.list \
    -O=${sample}.HC.vcf.gz
echo "gatk HaplotypeCaller generating VCF finished at $(date)"

# SNP mode
echo "SNP mode: gatk VariantRecalibrator started at $(date)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" VariantRecalibrator \
    -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
    -V ${sample}.HC.vcf.gz \
    --tmp-dir tmp \
    --resource:hap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=true,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=false,training=true,truth=true,prior=2.0 $GATK_bundle/dbsnp_146.hg38.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    --tranches-file ${sample}.HC.snps.tranches \
    --rscript-file ${sample}.HC.snps.plots.R \
    -O ${sample}.HC.snps.recal


gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" ApplyVQSR \
    -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
    -V ${sample}.HC.vcf.gz \
    --tmp-dir tmp \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ${sample}.HC.snps.tranches \
    --recal-file ${sample}.HC.snps.recal \
    -mode SNP \
    -O  ${sample}.HC.snps.VQSR.vcf.gz
echo "SNP mode: gatk VariantRecalibrator finished at $(date)"

# Indel mode
echo "Indel mode: gatk VariantRecalibrator started at $(date)"
gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" VariantRecalibrator \
    -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
    -V ${sample}.HC.snps.VQSR.vcf.gz \
    --tmp-dir tmp \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 ${GATK_bundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:dbsnp,known=false,training=true,truth=true,prior=2.0 ${GATK_bundle}/dbsnp_146.hg38.vcf.gz \
    -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL --max-gaussians 6 \
    --rscript-file ${sample}.HC.snps.indels.plots.R \
    --tranches-file ${sample}.HC.snps.indels.tranches \
    -O ${sample}.HC.snps.indels.recal


gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=15" ApplyVQSR \
    -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
    -V ${sample}.HC.snps.VQSR.vcf.gz \
    --tmp-dir tmp \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ${sample}.HC.snps.indels.tranches \
    --recal-file ${sample}.HC.snps.indels.recal \
    -mode INDEL \
    -O  VCF/${sample}.HC.VSQR.vcf.gz 
echo "Indel mode: gatk VariantRecalibrator finished at $(date)"
echo "GATK finished variant call and VCF files are stored at /work_dir/result/gatk/VCF"
rm -f *plots.R *recal*  *tranches *VQSR* *txt
rm ${sample}.*chr* ${sample}.input_variant_files.list