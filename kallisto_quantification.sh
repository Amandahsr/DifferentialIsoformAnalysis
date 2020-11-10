#Pipeline for running read quantification on kallisto.

#Initialise paired-end fastq files.
#4 files initialised because each read is split into 2.
fastq_read1_file1=read1_file1.fastq.gz
fastq_read1_file2=read1_file2.fastq.gz
fastq_read2_file1=read2_file1.fastq.gz
fastq_read2_file2=read2_file2.fastq.gz

#Parameters for kallisto. 
#Number of threads specified depends on number of files running.
threads=4
memory=80G
index=/mnt/projects/wlwtan/cardiac_epigenetics/pipeline/kallisto/kallisto_linux-v0.43.0/index/hg38/hg38
output_directory=directory_name

#Run kallisto. Queue as a binary job. 
#Bootstrapping is only neccessary if you would like to run packages that require bootstrapping (eg Sleuth).
#Queue different samples seperately. All files specified in same command is automatically taken as reads from the same sample. 
kallisto quant -t $threads -i $index -o $output_directory $fastq_read1_file1 $fastq_read1_file2 $fastq_read2_file1 $fastq_read2_file2  

