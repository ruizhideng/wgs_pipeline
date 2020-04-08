# command-line arguments, my work_dir = /home/ruizhi/data/wgs_pipeline/
# fastq raw data is stored at ${work_dir}/fasta_file
# assume your fastq filename format is: SAMPLE-ID_*.1.fastq.gz   SAMPLE-ID_*.2.fastq.gz
work_dir=$1

#reference
reference=${work_dir}/reference/genome/hg38.fa
GATK_bundle=${work_dir}/resources/gatk/bundle/hg38
# create result folders
if [ ! -d ${work_dir}/result/qc ]
then mkdir -p ${work_dir}/result/qc
fi

if [ ! -d ${work_dir}/result/cleanfq ]
then mkdir -p ${work_dir}/result/cleanfq
fi

if [ ! -d ${work_dir}/result/bwa ]
then mkdir -p ${work_dir}/result/bwa
fi

if [ ! -d ${work_dir}/result/gatk ]
then mkdir -p ${work_dir}/result/gatk
fi

if [ ! -d ${work_dir}/result/gatk/VCF ]
then mkdir -p ${work_dir}/result/gatk/VCF
fi

# fastqc: check fastq files quality, result stored at result/qc 
echo "fastqc checking data quality started at $(date)" 
cd ${work_dir}/result/qc
##fastqc ${work_dir}/fasta_file/*gz -o ${work_dir}/result/qc
echo "fastqc checking data quality finished at $(date)" 

# filter low_quality data and remove adapter
cd ${work_dir}/result/cleanfq
echo "trim_galore cutting adaptors started at $(date)" 
##ls ${work_dir}/fasta_file/*1.*.gz > ./1;cat 1
##ls ${work_dir}/fasta_file/*2.*.gz > ./2;cat 2
##paste 1 2 > config && rm 1 2
##cat ${work_dir}/result/cleanfq/config | while read id; do 
##    arr=(${id}); fq1=${arr[0]}; fq2=${arr[1]}; trim_galore -q 20 --phred33 --length 36 -e 0.1 --stringency 3 --paired $fq1 $fq2 -o ${work_dir}/result/cleanfq
##done
echo "trim_galore cutting adapters finished at $(date)" && rm config

echo "fastqc checking clean data quality started at $(date)" 
cd ${work_dir}/result/cleanfq
##fastqc ${work_dir}/result/cleanfq/*fq.gz -o ${work_dir}/result/qc
echo "fastqc checking clean data quality finished at $(date)"

# bwa mem: mapping reads to GRCh38
echo "bwa mem mapping reads started at $(date)"
cd ${work_dir}/result/bwa
ls ${work_dir}/result/cleanfq/*val_1.fq.gz > ./1;cat 1
ls ${work_dir}/result/cleanfq/*val_2.fq.gz > ./2;cat 2
paste 1 2 > config && rm 1 2
cat ${work_dir}/result/bwa/config | while read id; do
    arr=(${id}); fq1=${arr[0]}; fq2=${arr[1]}; tmp=`basename ${arr[0]}`; sample=`echo ${tmp%%_*}`; bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $reference $fq1 $fq2 | samtools view -bS - > $sample.bam
done
echo "bwa mapping reads finished at $(date)" && rm config

# sort bam files
echo "samtools sort sorting bam files started at $(date)"
ls *.bam | while read id; do samtools sort -@ 4 -m 500M -O bam -o ${work_dir}/result/bwa/`basename $id .bam`.sorted.bam $id;done
echo "samtools sort sorting bam files finished at $(date)"

#mark duplicated sequences
echo "gatk MarkDuplicates marks duplicated sequences started at $(date)"
cd ${work_dir}/result/bwa
ls *.sorted.bam| while read id; do
    gatk MarkDuplicates \
        -I $id \
        --REMOVE_DUPLICATES=true \
        --VALIDATION_STRINGENCY=SILENT \
        -M ${work_dir}/result/gatk/`basename $id .sorted.bam`.markup_metrics.txt \
        -O ${work_dir}/result/gatk/`basename $id .sorted.bam`.sorted.markdup.bam 
done     
echo "gatk MarkDuplicates marks duplicated sequences finished at $(date)"

# samtools index building markdup.bam index 
cd ${work_dir}/result/gatk/
echo "samtools index building markdup.bam index started at $(date)"
ls *markdup.bam| while read id; do samtools index $id; done
echo "samtools index building markdup.bam index finished at $(date)"
