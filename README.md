# This pipeline used to call variants, no biological duplication.
# Please strictly follow my steps and samples naming
# 1. Command-line argument is work_dir

1. work_dir=$1  # $1 that meas your command-line argument, my work_dir is: /home/ruizhi/data/wgs_pipeline/

   reference=/home/ruizhi/data/wgs_pipeline/reference/genome # hg38 is stored here

   GATK_bundle=/home/ruizhi/data/wgs_pipeline/resources/gatk/bundle/hg38 # GATK bundle is stored here

2. your fasta files are stored at ${work_dir}/fasta_file, like /home/ruizhi/data/wgs_pipeline//fasta_file

   And please put all of your fasta files to ${work_dir}/fasta_file
   
   We will use bwa to map reads to hg38.

3. GATK bundle is stored at ${GATK_bundle}, like /home/ruizhi/data/wgs_pipeline/resources/gatk/bundle/hg38
   
   And please download all of your GATK bundle to ${GATK_bundle}
   
   GATK bundle include all resource that we need for gatk variant calling.
   
   You can find instruction for how to download those data. Please follow my instruction to download the same GATK bundle, 
   
   then you don't need to build any index for those gatk bundle.

4. This script just used to analyze non-biological samples. Your fasta files are pair-end sequecing.

   The format of sample names must be: SAMPLE-ID_*.1.fq.gz      SAMPLE-ID_*.2.fq.gz
   
   like KPGP-00246_L001_R1.fq.gz       KPGP-00246_L001_R2.fq.gz
   
        KPGP-00246_*1.fq.gz             KPGP-00246_*2.fq.gz
        
   Please strictly follow my naming rule!!!

# 2. Preparation

1. Please USING conda to intall softwares:

   conda install trim_galore
   
   conda install bwa
   
   conda install samtools   
 
   conda install gatk/4.1.2.0 (if you can't conda gatk/4.1.2.0, please intall by source, and add to enciroment varible)

2. Please download human reference genome hg38 to reference folder, like /home/ruizhi/data/wgs_pipeline/reference/genome

   cd ${reference} 
  
   wget http://hgdownload.cse.ucsc.edu/goldenPath.hg38/bigZips/hg38.fa.gz

3. please download GATK bundle and stored at ${GATK_bundle}. (recommend use FileZilla)

   GATK bundle can be downloaded from ftp-site (LEAVE THE PASSWORD BLANK):

   ftp.broadinstitute.org/bundle

   username: gsapubftp-anonymous

   password:
   
# 3. RUN script

1. Firstly, run preparation.sh

   ./preparation.sh /work_dir/
   
2. Then, qc and mapping, get bam files

   ./fastq_to_bam.sh /work_dir/

3. Last, call variant
   
   ./GATK_variant_call.sh /work_dir/
   


