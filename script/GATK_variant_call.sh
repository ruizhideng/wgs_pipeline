# command-line arguments, my work_dir = /home/ruizhi/data/wgs_pipeline/
# assume your sorted bam filename format is: SAMPLE-ID.sorted.bam
# assum your sorted bam index filename format is: SAMPLE-ID.sorted.markup.bam
work_dir=$1

#reference
reference=${work_dir}/reference/genome/hg38.fa
GATK_bundle=${work_dir}/resources/gatk/bundle/hg38

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

# excute BQSR
cd ${work_dir}/result/gatk/
echo "gatk BaseRecalibrator recalibrating started at $(date)"
ls *.sorted.markdup.bam| while read id; do
    gatk BaseRecalibrator \
        -R $GATK_bundle/Homo_sapiens_assembly38.fasta \
        -I ${work_dir}/result/gatk/$id \
        --known-sites $GATK_bundle/dbsnp_146.hg38.vcf.gz \
        --known-sites $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        -O `basename $id .bam`.recal_data.table
done
echo "gatk BaseRecalibrator recalibrating finished at $(date)"

# gatk ApplyBQSR 
echo "gatk ApplyBQSR started at $(date)"
ls *sorted.markdup.bam| while read id; do
    gatk ApplyBQSR \
        --bqsr-recal-file `basename $id .bam`.recal_data.table \
        -R $GATK_bundle/Homo_sapiens_assembly38.fasta \
        -I $id \
        -O `basename $id .bam`.BQSR.bam
done
echo "gatk ApplyBQSR finished at $(date)"

# build *sorted.markdup.BQSR.bam index
echo "samtools index building index started at $(date)"
ls *sorted.markdup.BQSR.bam| while read id; do samtools index $id; done
echo "samtools index building index finished at $(date)"

# generate raw VCF files
echo "gatk HaplotypeCaller generating VCF started at $(date)"
ls *sorted.markdup.BQSR.bam| while read id; do
    gatk HaplotypeCaller \
	    -R $GATK_bundle/Homo_sapiens_assembly38.fasta \
	    -I $id \
	    -O `basename $id .sorted.markdup.BQSR.bam`.HC.vcf.gz 
done
echo "gatk HaplotypeCaller generating VCF finished at $(date)"

# SNP mode
echo "SNP mode: gatk VariantRecalibrator started at $(date)"
ls *.HC.vcf.gz| while read id; do
    gatk VariantRecalibrator \
        -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
        -V ${id} \
        --resource:hap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
        --resource:1000G,known=false,training=true,truth=true,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        --resource:dbsnp,known=false,training=true,truth=true,prior=2.0 $GATK_bundle/dbsnp_146.hg38.vcf.gz \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
        --tranches-file `basename $id .HC.vcf.gz`.HC.snps.tranches \
        --rscript-file `basename $id .HC.vcf.gz`.HC.snps.plots.R \
        -O `basename $id .HC.vcf.gz`.HC.snps.recal
done

ls *.HC.vcf.gz| while read id; do
    gatk ApplyVQSR \
        -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
        -V $id \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file `basename $id .HC.vcf.gz`.HC.snps.tranches \
        --recal-file `basename $id .HC.vcf.gz`.HC.snps.recal \
        -mode SNP \
        -O  `basename $id .HC.vcf.gz`.HC.snps.VQSR.vcf.gz
done
echo "SNP mode: gatk VariantRecalibrator finished at $(date)"

# Indel mode
echo "Indel mode: gatk VariantRecalibrator started at $(date)"
ls *.HC.snps.VQSR.vcf.gz| while read id; do
    gatk VariantRecalibrator \
        -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
        -V $id \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 ${GATK_bundle}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --resource:dbsnp,known=false,training=true,truth=true,prior=2.0 ${GATK_bundle}/dbsnp_146.hg38.vcf.gz \
        -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode INDEL \
        --rscript-file `basename $id .HC.snps.VQSR.vcf.gz`.HC.snps.indels.plots.R \
        --tranches-file `basename $id .HC.snps.VQSR.vcf.gz`.HC.snps.indels.tranches \
        -O `basename $id .HC.snps.VQSR.vcf.gz`.HC.snps.indels.recal
done

ls *.HC.snps.VQSR.vcf.gz| while read id; do
    gatk ApplyVQSR \
        -R ${GATK_bundle}/Homo_sapiens_assembly38.fasta \
        -V $id \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file `basename $id .HC.snps.VQSR.vcf.gz`.HC.snps.indels.tranches \
        --recal-file `basename $id .HC.snps.VQSR.vcf.gz`.HC.snps.indels.recal \
        -mode INDEL \
        -O  ${work_dir}/result/gatk/VCF/`basename $id .HC.snps.VQSR.vcf.gz`.HC.VSQR.vcf.gz 
done
echo "Indel mode: gatk VariantRecalibrator finished at $(date)"
echo "GATK finished variant call and VCF files are stored at /work_dir/result/gatk/VCF"
rm -f *plots.R *recal*  *tranches *VQSR* *txt *.sorted.makdup.bam*