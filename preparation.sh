#!/usr/bin/bash

# my work_dir=/home/ruizhi/data/wgs_pipeline
# fastq raw data is stored at ${work_dir}/test_data
work_dir=$1
cd ${work_dir}

if [ ! -d ${work_dir}/test_data ]
then mkdir -p ${work_dir}/test_data
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

#download test data
cd ${work_dir}/test_data
nohup wget -c -r -nd -np -k -L -p ftp://biodisk.org/Release/KPGP/KPGP_Data_2015_Release_Candidate/WES/KPGP-00245/ 1>/dev/null 2>&1 &

# download reference GRCh38
cd ${work_dir}/reference/genome
nohup wget http://hgdownload.cse.ucsc.edu/goldenPath.hg38/bigZips/hg38.fa.gz &
gzip hg38.fa.gz

# download gatk bundle hg38
# ftp://ftp.broadinstitute.org/bundle/hg38/
cd ${work_dir}/resource/gatk/bundle
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz &
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai & nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi &


# build reference index
cd ${work_dir}/index/bwa
nohup time bwa index -a bwtsw -p /home/ruizhi/data/wgs_pipeline/reference/genome/hg38.fa /home/ruizhi/data/wgs_pipeline/reference/genome/hg38.fa  1>hg38.bwa_index.log 2>&1 &
  

