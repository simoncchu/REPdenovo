__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
import shutil
import glob
from subprocess import *


#calculate coverage 
def calcCoverage(fpath, read_length, genome_length, VERBOSE):
    #print fpath
    cnt_lines=0
    cmd=""
    if fpath.lower().endswith(('.fq', '.fastq')):
        cmd="echo $(wc -l {0})".format(fpath)
    elif fpath.lower().endswith(('.fastq.gz', 'fq.gz', 'gz')):
        cmd="echo $(zcat {0} | wc -l)".format(fpath)
        #p=Popen(cmd, shell = True, stdout = PIPE)
        #p.communicate()
        #cnt_lines=int(p.stdout.readline())
    else:
        print "Something wrong with the raw reads files format:", fpath

    if VERBOSE != 0:
        print "Running command: "+ cmd +"..."
    tp=tuple(Popen(cmd, shell = True, stdout = PIPE).communicate())
    lcnt=str(tp[0]).split()
    print "The number of reads in file {0} is {1}".format(fpath,tp)
    cnt_lines=int(lcnt[0])
    cnt_reads=int(cnt_lines)/4
    cov=float(cnt_reads*read_length)/float(genome_length)
    return cov

#calcCoverage(sys.argv[1], 0, 0)

#concatenate files of given suffix
#Input: suffix
#Output:
'''
def concatenateFiles(suffix, sout):
    with open(sout, 'wb') as outfile:
        for filename in glob.glob(suffix):
            with open(filename, 'rb') as readfile:
                shutil.copyfileobj(readfile, outfile)

'''


'''
Description:
    Read in a fasta file
Parameters:
    sffa: path of fa file
    brstrip: whether strip '\n' or others can be striped by rstrip()
'''
def readContigFa(sffa,brstrip):
    dcontigs={}
    name=""
    seq=""
    last_name=""
    ffa=open(sffa)
    for line in ffa:
        if brstrip==True:
            line=line.rstrip()
        if line[0]=='>':
            name=line[1:]
            if last_name!="":
                dcontigs[last_name]=seq
                seq=""
            last_name=name
        else:
            seq=seq+line
    dcontigs[last_name]=seq
    ffa.close()
    return dcontigs

