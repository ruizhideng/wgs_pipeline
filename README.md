# This pipeline used to call variants
# 1. Command-line argument is work_dir

1. work_dir is path/, like /your_path/wgs_analysis

2. your fasta files are stored at work_dir/fasta_files

3. GATK bundle is stored at /work_dir/resources/gatk/bundle

# 2. Preparation

1. Please USING conda to intall softwares:

   trim_galore, bwa, samtools, gatk/4.0

2. Please download human reference genome hg38.

   wget http://hgdownload.cse.ucsc.edu/goldenPath.hg38/bigZips/hg38.fa.gz

3. please download GATK bundle and stored at /work_dir/resources/gatk/bundle. (recommend use FileZilla)

   GATK bundle can be downloaded from ftp-site (LEAVE THE PASSWORD BLANK):

   ftp.broadinstitute.org/bundle

   username: gsapubftp-anonymous

   password:
   
# 3. RUN script

1. Firstly, run preparation.sh

   ./preparation.sh /work_dir/
   
2. Then, qc and mapping, get bam files
   ./


