#!/usr/bin/bash

# my work_dir=/home/ruizhi/data/wgs_pipeline
# fastq raw data is stored at ${work_dir}/fasta_file
work_dir=$1
reference=${work_dir}/reference/genome/hg38.fa

if [ ! -d ${work_dir}/result/platypus ]
then mkdir -p ${work_dir}/result/platypus
fi 

# Please download Platypus from https://github.com/RahmanTeamDevelopment/Platypus
# install Platypus instruction, before install Platupus, you have to first install Virttualenv and Python (version 2.7.X is currently supported)
conda install Virttualenv
# cd /Platypus
# ./install.sh
##### Be careful! when you install platypus, your python will be downgraded to 2.7 version!!!
env/bin/platypus callVariants --bamFiles=${work_dir}/result/bwa/ 
--refFile=reference.fa \
--output=variant_call.vcf \
--logFileName=platypus_log_file.txt \
--minPosterior=20 \
