# my work_dir=/home/ruizhi/data/wgs_pipeline
# my platypus=/home/ruizhi/data/biosoft/Platypus-master/env/bin/platypus
# bam files are stored at ${work_dir}/result/bwa
work_dir=$1
platypus=$2
platypus_reference=${work_dir}/resources/gatk/bundle/hg38/Homo_sapiens_assembly38.fasta

if [ ! -d ${work_dir}/result/platypus ]
then mkdir -p ${work_dir}/result/platypus
fi 

# Please download Platypus from https://github.com/RahmanTeamDevelopment/Platypus
# install Platypus instruction, before install Platupus, you have to first install Virttualenv and Python (version 2.7.X is currently supported)
# conda install Virttualenv
# cd /Platypus
# ./install.sh
##### Be careful! when you install platypus, your python will be downgraded to 2.7 version!!!
cd ${work_dir}/result/gatk/
echo "platypus started at $(date)"
ls *markdup.bam | while read id; do
    ${platypus} callVariants --bamFiles=$id --refFile=${platypus_reference} --output=${work_dir}/result/platypus/`basename $id .markdup.bam`.vcf --logFileName=`basename $id .sorted.bam`.platypus_log_file.txt --minPosterior=5
done
echo "platypus finished at $(date)"