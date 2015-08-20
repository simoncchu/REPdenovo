__author__ = 'Chong Chu'

import sys
import os
from subprocess import *

#sam reads sorted by read name
def filterPEByMapQuality(sfsam,mapq, sfsam_output):
    print "Filtering out PE reads with mapping quality are smaller than "+str(mapq)+" ..."

    #first check whether files are existing
    if os.path.exists(sfsam)==False:
        print sfsam+" is not existing."
        return

    fin_reads=open(sfsam,'r')
    fout_reads=open(sfsam_output,'w')

    pre_id=""
    cur_id=""
    pre_mapq=0
    cur_mapq=0
    is_first=True
    pre_line=""

    for line in fin_reads:
        #bypass head
        if line[0]=="@":
            fout_reads.write(line)
            continue

        parts=line.split()
        if(len(parts)<11):
            print "Wrong reads!!!"
            continue
        else:
            cur_id=parts[0] #qname
            cur_mapq=int(parts[4]) #mapq
            if is_first==True:
                pre_id=cur_id
                pre_mapq=cur_mapq
                is_first=False
                pre_line=line
                continue
            else:
                #check whether is a pair and also the mapping quality of the pair is larger than threshold
                if (cur_id==pre_id) and (cur_mapq>=mapq) and (pre_mapq>=mapq):
                    #write into file
                    fout_reads.write(pre_line)
                    fout_reads.write(line)

        pre_id=cur_id
        pre_mapq=cur_mapq
        pre_line=line

    fout_reads.close()
    fin_reads.close()

#def filterReadsByMapQuality(sfsam,mapq, sfsam_output):


##filter by cigar
def filterByCigar():
    return


def filterSam(sfsam, mapq,sfsam_output):
    print "First, filter out those unmapped reads in {0} ...".format(sfsam)
    #First, only keep those fully mapped reads
    sintermediate="{0}_intermediate.sam".format(sfsam)
    cmd="samtools view -h -S -F 4 {0} > {1}".format(sfsam,sintermediate)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools view -h -S -b {0} > {1}.bam".format(sintermediate,sintermediate)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    os.remove(sfsam)

    ##keep the mapped reads(for calculating coverage), covert to bam and sort
    sbam_for_cov="{0}_for_coverage".format(sfsam)
    print "Sort filtered {0} ...".format(sintermediate)
    cmd="samtools sort {0}.bam {1}.sorted".format(sintermediate,sbam_for_cov)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools index {0}.sorted.bam".format(sbam_for_cov)

    print "Sort filtered {0} by read name ...".format(sfsam)
    #first covert to bam, then sort
    cmd="samtools sort -n {0}.bam {1}.sortbyname".format(sintermediate,sintermediate)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    # then covert back to sam
    cmd="samtools view -h {0}.sortbyname.bam > {1}".format(sintermediate,sintermediate)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    print "Filter out pairs that mapping quality are low, or only one in pair is qualified...."
    #Then, filter through mapping quality
    filterPEByMapQuality(sintermediate, mapq, sfsam_output)

    #remove useless files
    os.remove(sintermediate)
    #os.remove("{0}.bam".format(sintermediate))
    os.remove("{0}.sortbyname.bam".format(sintermediate))
    #os.remove("{0}.sortbyname.bam.bai".format(sintermediate)) 
