# my work_dir=/home/ruizhi/data/wgs_pipeline
# fastq raw data is stored at ${work_dir}/fasta_file
work_dir=$1
genome=${work_dir}/reference/genome
GATK_bundle=${work_dir}/resources/gatk/bundle/hg38
cd ${work_dir}

if [ ! -d ${work_dir}/fasta_file ]
then mkdir -p ${work_dir}/fasta_file
fi

if [ ! -d ${work_dir}/reference/genome ]
then mkdir -p ${work_dir}/reference/genome 
fi

if [ ! -d ${work_dir}/index/bwa ]
then mkdir -p ${work_dir}/index/bwa  
fi

if [ ! -d ${work_dir}/resources/gatk/bundle ]
then mkdir -p ${work_dir}/resources/gatk/bundle 
fi

# download test data
##cd ${work_dir}/fasta_file
##nohup wget -c -r -nd -np -k -L -p ftp://biodisk.org/Release/KPGP/KPGP_Data_2015_Release_Candidate/WES/KPGP-00245/ 1>/dev/null 2>&1 &

# download reference GRCh38 and build reference index
cd $genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz
echo "bwa index hg38 started at $(date)"
bwa index -a bwtsw $genome/hg38.fa
echo "bwa index hg38 finished at $(date)" 
# samtools index building hg38 index
echo "samtools index hg38 started at $(date)"
samtools faidx $genome/hg38.fa
echo "samtools index hg38 finished at $(date)"  
# download gatk bundle hg38 from ftp-sites
# ftp://ftp.broadinstitute.org/bundle/hg38/
# username: gsapubftp-anonymous
# password:
# build GATK_bundle index
cd $GATK_bundle
gunzip $GATK_bundle/Homo_sapiens_assembly38.fasta.gz