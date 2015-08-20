__author__ = 'Chong Chu'

import sys
import os
from subprocess import *

'''
sh removeDupRepeatsOfOneContigSet.sh $1 $1.no_dup.fa $2
sh mergeContigs.sh $1.no_dup.fa
sh removeDupRepeatsOfOneContigSet.sh $1.no_dup.fa.merged.fa $1.no_dup.fa.merged.fa.no_dup.fa $3
sh removeContainedRepeatsOfOneContigSet.sh $1.no_dup.fa.merged.fa.no_dup.fa $1.no_dup.fa.merged.fa.no_dup.no_contain.fa $3
sh alignContigs2Bcmk.sh $1.no_dup.fa.merged.fa.no_dup.fa
'''

def removeDuplicateContained(fcontig, foutput, cutoff, brm_contained):
    #remove duplicate or contained contigs
    cmd="samtools faidx {0}".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="./TERefiner_1 -U -r {0} -o {1}".format(fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="samtools faidx {0}".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="bwa index {0}".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="bwa mem -a {0} {1} > {2}.itself.sam".format(fcontig,fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    cmd="samtools view -h -S -b {0}.itself.sam > {1}.itself.bam".format(fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools sort {0}.itself.bam {1}.itself.sort".format(fcontig,fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools index {0}.itself.sort.bam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    if brm_contained==True:
        cmd="./TERefiner_1 -P -b {0}.itself.sort.bam -r {1} -o {2} -c {3} -g".format(fcontig,fcontig,foutput,cutoff)
    else:
        cmd="./TERefiner_1 -P -b {0}.itself.sort.bam -r {1} -o {2} -c {3}".format(fcontig,fcontig,foutput,cutoff)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    ##clean all the temporary files
    cmd="rm {0}.sa".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.pac".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.bwt".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.ann".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.amb".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.sam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.bam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.sort.bam".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.itself.sort.bam.bai".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="rm {0}.fai".format(fcontig)
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def mergeContigs(fout_folder, cutoff_dup_bf_merge, cutoff_dup_af_merge):
    fcontig=fout_folder+"contigs.fa"
    foutput=fcontig+"_no_dup.fa"
    removeDuplicateContained(fcontig, foutput, cutoff_dup_bf_merge, False)

    cmd="./ContigsMerger -s 0.2 -i1 -6.0 -i2 -6.0 -x 15 -k 10 -o {0}.merge.info {1} > {2}.merged.fa".format(foutput,foutput,foutput)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    fcontig="{0}.merged.fa".format(foutput)
    foutput="{0}.no_dup.fa".format(fcontig)
    removeDuplicateContained(fcontig,foutput, cutoff_dup_af_merge ,False)
    fcontig=foutput
    foutput=fcontig+".no_contained.fa"
    removeDuplicateContained(fcontig,foutput, cutoff_dup_af_merge ,True)

    #rename contigs.fa
    os.rename(fout_folder+"contigs.fa", fout_folder+"original_contigs_before_merging.fa")
    os.rename(foutput,fout_folder+"contigs.fa")

