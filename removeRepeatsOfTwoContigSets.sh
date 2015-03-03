#Given two files XX1.fa and XX2.fa, align to each other
#Step 0. Remove duplicate contigs
samtools faidx $1
./TERefiner_1 -U -r $1 -o $1
samtools faidx $2
./TERefiner_1 -U -r $2 -o $2
#
#Step 1. Index both XX1.fa and XX2.fa
bwa index $1
bwa index $2
samtools faidx $1
samtools faidx $2
#Step 2. Convert to fastq files
awk -f fasta2fastq.awk $1 > $1.fq
awk -f fasta2fastq.awk $2 > $2.fq
#
#Step 3. align XX1.fq to XX2.fa, XX2.fq to XX1.fa
bwa mem $2 $1.fq > $1.sam
bwa mem $1 $2.fq > $2.sam
#
#Step 4. convert to bam, index and sort
samtools view -h -S -b $1.sam > $1.bam
samtools view -h -S -b $2.sam > $2.bam
samtools sort $1.bam $1.sort
samtools sort $2.bam $2.sort
samtools index $1.sort.bam
samtools index $2.sort.bam
#
#Step 5. Find out repeat ones 
echo "Repeat contigs that need to remove in "$1":"
./TERefiner_1 -T -b $1.sort.bam -r $2 -s $1 -o $3 -c $5
echo "--------------------------------------------"
echo "Repeat contigs that need to remove in "$2":"
./TERefiner_1 -T -b $2.sort.bam -r $1 -s $2 -o $4 -c $5