__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
from subprocess import *

#get file size
def getSize(filename):
    st = os.stat(filename)
    return st.st_size

#####################
#output:
# -1, if no output
#0, otherwise
####################
def collectHighFreqKmers(kmer_freq, fkmer, fname_out):
    sss="Collect high frequency (larger than {0}) kmers...".format(kmer_freq)
    print sss
    cnt=0
    index=0
    '''
    with open(fkmer, 'rt') as fin:
        with open(fname_out, 'wt') as fout:
            for stemp in fin:
                if stemp[0]=='>':
                    cnt=int(stemp[1:])
                elif cnt>=kmer_freq:
                    read_name="@seq{0}\n".format(index)
                    index+=1
                    fout.write(read_name)
                    fout.write(stemp)
                    fout.write("+\n")
                    seq_len=len(stemp)
                    iss=1
                    quality=""
                    while iss<seq_len:
                        quality+="5"
                        iss+=1
                    quality+="\n"
                    fout.write(quality)
    '''
    fin=open(fkmer, 'rt')
    fout=open(fname_out, 'wt')
    for stemp in fin:
        if stemp[0]=='>':
            cnt=int(stemp[1:])
        elif cnt>=kmer_freq:
            read_name="@seq{0}\n".format(index)
            index+=1
            fout.write(read_name)
            fout.write(stemp)
            fout.write("+\n")
            seq_len=len(stemp)
            iss=1
            quality=""
            while iss<seq_len:
                quality+="5"
                iss+=1
            quality+="\n"
            fout.write(quality)
    fout.close()
    fin.close()

    if getSize(fname_out)<=0:
        ss="No kmer found when frequency is {0}".format(kmer_freq)
        print ss
        return -1
    return 0


#####################
#output:
# -1, if no output
#0, otherwise
####################
def collectHighFreqKmersOfRegion(kmer_freq, slack ,fkmer, fname_out):
    sss="Collect high frequency (larger than {0}) kmers...".format(kmer_freq)
    print sss
    cnt=0
    index=0
    fin=open(fkmer, 'rt')
    fout=open(fname_out, 'wt')
    for stemp in fin:
        if stemp[0]=='>':
            cnt=int(stemp[1:])
        elif cnt>=kmer_freq:
            read_name="@seq{0}\n".format(index)
            index+=1
            fout.write(read_name)
            fout.write(stemp)
            fout.write("+\n")
            seq_len=len(stemp)
            iss=1
            quality=""
            while iss<seq_len:
                quality+="5"
                iss+=1
            quality+="\n"
            fout.write(quality)
    fout.close()
    fin.close()

    if getSize(fname_out)<=0:
        ss="No kmer found when frequency is {0}".format(kmer_freq)
        print ss
        return -1
    return 1


#Input:
    #k_len: usually is the kmer length minus one.
def assembleKmer(vpath, k_len, ffastq, min_contig_len, output_folder, VERBOSE):
    print "Assembly high frequency kmers ..."
    if vpath!="" and vpath[-1]!="/":
        vpath+="/"
    vpath_h = vpath+"velveth"
    vpath_g = vpath+"velvetg"
    #print vpath_h########################################################################################################
    assert os.path.exists(vpath_h),"velveth is not installed, or given the wrong path in configuration file!!!"
    assert os.path.exists(vpath_g),"velvetg is not installed, or given the wrong path in configuration file!!!"

    newpath=output_folder
    #print newpath######################################################################################################################################
    if not os.path.exists(newpath): os.makedirs(newpath)

    cmd="{0} {1} {2} -fastq -short {3}".format(vpath_h,newpath,k_len,ffastq)
    if VERBOSE != 0:
        print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} {1} -min_contig_lgth {2} -unused_reads yes".format(vpath_g,newpath,min_contig_len)
    if VERBOSE != 0:
        print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()

#/scratch2/chongchu/velvet-master/velveth ./temp_assembly 49 -fastq  -short kmer_NA12878_k50_morethan5k.fq
#/scratch2/chongchu/velvet-master/velvetg temp_NA12878_k50_5k/ -min_contig_lgth 100

'''
def removeDuplicate(fa_old, fa_new):
    print "Generate the unique fa..."

    with open(fa_old, 'rt') as fin:
        with open(fa_new, 'wt') as fout:

'''

'''
#Given one files XX1.fa , align to itself
#Step 0. Remove exact duplicate ones
samtools faidx $1
TERefiner_1 -U -r $1 -o $1
#Step 1. Index XX1.fa
bwa index $1
samtools faidx $1
#Step 2. Convert to fastq files
awk -f fasta2fastq.awk $1 > $1.fq

#Step 3. align XX1.fq to XX1.fa
bwa mem -a -M $1 $1.fq > $1.itself.sam

#Step 4. convert to bam, index and sort
samtools view -h -S -b $1.itself.sam > $1.itself.bam
samtools sort $1.itself.bam $1.itself.sort
samtools index $1.itself.sort.bam

#Step5. Find out repeat ones
TERefiner_1 -O -b $1.itself.sort.bam -r $1 -o $2 -l $3
'''
#def rmDupOneContigSet():
