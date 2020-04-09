# cnvator CNV
work_dir=$1
cnvnator=$2 #/PATH/TO/cnvnator/
reference=${work_dir}/reference/genome/hg38.fa

if [ ! -d ${work_dir}/result/CNV/CNVnator ]
then mkdir -p ${work_dir}/result/CNV/CNVnator
fi

# Extract read mapping
echo "CNVnator calling CNV started at $(date)"
cd ${work_dir}/result/gatk/
ls *sorted.markdup.bam| while read id; do
    $cnvnator -root ${work_dir}/result/CNV/CNVnator/`basename ${id} .sorted.markdup.bam`.root -tree ${id}
done

# Generate histogram
cd ${work_dir}/result/CNV/CNVnator
ls *.root| while read id; do
    $cnvnator -root ${id} -his 1000 -fasta ${reference}
done

# Calculate statistics
ls *.root| while read id; do 
    $cnvnator -root ${id} -stat 1000
done

#Partition
ls *.root| while read id; do 
    $cnvnator -root ${id} -partition 1000
done

# Call CNVs
ls *.root| while read id; do 
    $cnvnator -root ${id} -call 1000 > `basename ${id} .root`.call.txt
done
echo "CNVnator calling CNV finished at $(date)"