# Manta CNV tools
work_dir=$1
manta=$2 #/PATH/TO/configManta.py/ like /home/ruizhi/data/biosoft/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py 
reference=${work_dir}/reference/genome/hg38.fa
# create result folders
if [ ! -d ${work_dir}/result/CNV/Manta ]
then mkdir -p ${work_dir}/result/CNV/Manta
fi

# Manta calling
echo "Manta calling to CNVs started at $(date)"
cd ${work_dir}/result/gatk/
ls *sorted.markdup.bam| while read id; do
    $manta --bam $id \
    --referenceFasta ${reference} \
    --runDir ${work_dir}/result/CNV/Manta
done
echo "Manta calling to CNVs finished at $(date)"
