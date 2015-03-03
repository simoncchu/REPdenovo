#Given one files XX1.fa , align to itself
#Step 0. Remove exact duplicate ones
samtools faidx $1
./TERefiner_1 -U -r $1 -o $1
#Step 1. Index XX1.fa
bwa index $1
samtools faidx $1
#Step 2. Convert to fastq files
awk -f fasta2fastq.awk $1 > $1.fq
#
#Step 3. align XX1.fq to XX1.fa
bwa mem -a -M $1 $1.fq > $1.itself.sam
#
#Step 4. convert to bam, index and sort
samtools view -h -S -b $1.itself.sam > $1.itself.bam
samtools sort $1.itself.bam $1.itself.sort
samtools index $1.itself.sort.bam
#
#Step5. Find out repeat ones
./TERefiner_1 -O -b $1.itself.sort.bam -r $1 -o $2 -c $3