# GRIDSS CNV tools
# -t meanns threads
work_dir=$1
gridss=$2 #/PATH/TO/gridss/ like /home/ruizhi/data/biosoft/gridss/, gridss installed directory
reference=${work_dir}/reference/genome/hg38.fa

if [ ! -d ${work_dir}/result/CNV/gridss ]
then mkdir -p ${work_dir}/result/CNV/gridss 
fi

echo "gridss calling CNV started at $(date)"
cd ${work_dir}/result/gatk/
ls *sorted.markdup.bam| while read id; do
    ${gridss}/gridss.sh  \
    -t 1 \
    -r ${reference} \
    -o ${work_dir}/result/CNV/gridss/`basename ${id} .sorted.markdup.bam`.vcf.gz \
    -a ${work_dir}/result/CNV/gridss/`basename ${id} .sorted.markdup.bam`.assembly.bam \
    -j ${gridss}/gridss-2.8.3-gridss-jar-with-dependencies.jar \
    -w ${work_dir}/result/CNV/gridss \
    --jvmheap 4g \
    $id;done
echo "gridss calling CNV finished at $(date)"
