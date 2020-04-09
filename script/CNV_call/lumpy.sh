#lumpy CNV tool
work_dir=$1
scripts=$2 # /path/to/scripts, like /home/ruizhi/data/biosoft/lumpy-sv/scripts
           # stored with extractSplitReads_BwaMem, pairend_distro.py 
reference=${work_dir}/reference/genome/hg38.fa
# create result folders
if [ ! -d ${work_dir}/result/CNV/lumpy ]
then mkdir -p ${work_dir}/result/CNV/lumpy
fi

if [ ! -d ${work_dir}/result/CNV/lumpy/VCF ]
then mkdir -p ${work_dir}/result/CNV/lumpy/VCF
fi
# preparation for lumpy
# extract discordant paired-end alighments
echo "lumpy preparing started at $(date)"
cd ${work_dir}/result/gatk
ls *.sorted.markdup.bam| while read id; do
    samtools view -bh -F 1294 ${id} > ${work_dir}/result/CNV/lumpy/`basename ${id} .sorted.markdup.bam`.discordants.bam
done
#extract the split-read alighments
ls *.sorted.markdup.bam| while read id; do
    samtools view -h ${id} \
        | ${scripts}/extractSplitReads_BwaMem -i stdin \
        | samtools view -Sb - \
        > ${work_dir}/result/CNV/lumpy/`basename ${id} .sorted.markdup.bam`.splitters.bam
done
# sort both alighments
cd ${work_dir}/result/CNV/lumpy
ls *.discordants.bam| while read id; do
    samtools sort -@ 4 ${id} \
    -o `basename ${id} .discordants.bam`.discordants.sorted.bam
done
ls *.splitters.bam| while read id; do
    samtools sort -@ 4 ${id} \
    -o `basename ${id} .splitters.bam`.splitters.sorted.bam
done
echo "lumpy preparation finished at $(date)"

# lumpy CNVs
echo "lumpy started at $(date)"
# First, generate empirical insert size statistics on each library in the BAM file
cd ${work_dir}/result/gatk
ls *.sorted.markdup.bam| while read id; do
    samtools view -r readgroup1 ${id} \
    | tail -n+100000 \
    | ${scripts}/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 10000 \
    -o ${work_dir}/result/CNV/lumpy/`basename ${id} .sorted.markdup.bam`.lib1.histo
done

# run lympy with paried-end and split-reads
cd ${work_dir}/result/CNV/lumpy/
ls *.discordants.sorted.bam| while read id; do
    lumpy -mw 4 -tt 0 \
    -pe id:`basename ${id} .discordants.sorted.bam`,bam_file:${id},histo_file:`basename ${id} .discordants.sorted.bam`.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:`basename ${id} .discordants.sorted.bam`,bam_file:`basename ${id} .discordants.sorted.bam`.splitters.sorted.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    > ${work_dir}/result/CNV/lumpy/VCF/`basename ${id} .discordants.sorted.bam`.vcf
done
echo "lumpy finished at $(date)"