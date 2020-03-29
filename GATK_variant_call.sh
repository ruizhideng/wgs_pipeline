# command-line arguments, my work_dir = /home/ruizhi/data/wgs_pipeline/
# assume your sorted bam filename format is: SAMPLE-ID.sorted.bam
# assum your sorted bam index filename format is: SAMPLE-ID.sorted.markup.bam
work_dir=$1

#reference
reference=/home/ruizhi/data/wgs_pipeline/reference/genome/hg38.fa
GATK_bundle=/home/ruizhi/data/wgs_pipeline/resources/gatk/bundle

#mark duplicated sequences
echo "gatk MarkDuplicates marks duplicated sequences started at $(date)"
cd ${work_dir}/result/bwa
ls *.sorted.bam| while read id; do
    gatk MarkDuplicates \
        -I $id \
        -M `basename $id .sorted.bam`.markup_metrics.txt \
        -O ${work_dir}/result/gatk/`basename $id .sorted.bam`.sorted.markdup.bam 
done     
echo "gatk MarkDuplicates marks duplicated sequences finished at $(date)"

# samtools index building markdup.bam index 
##cd ${work_dir}/result/gatk/
##echo "samtools index building markdup.bam index started at $(date)"
##ls *markdup.bam| while read id; do samtools index $id; done
##echo "samtools index building markdup.bam index finished at $(date)"


