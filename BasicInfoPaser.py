__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
import shutil
import glob
from subprocess import *
from Utility import OUTPUT_FOLDER
from Utility import get_output_folder
from Utility import print_command

def calc_total_coverage(bpaired,sfleft_reads,sfright_reads,sfsingle_reads,READ_LENGTH,GENOME_LENGTH):
    print "Calculating average coverage...."
    OUTPUT_FOLDER=get_output_folder()
    sfcoverage=OUTPUT_FOLDER+"reads_coverage.txt"
    f=2.0

    if os.path.exists(sfcoverage)!=True:
        #print "Test: ", sfleft_reads###########################################################################################################3
        if bpaired==True:
            f=calc_coverage(sfleft_reads,READ_LENGTH, GENOME_LENGTH) + calc_coverage(sfright_reads,READ_LENGTH, GENOME_LENGTH)
        else:
            f=calc_coverage(sfsingle_reads,READ_LENGTH, GENOME_LENGTH)
        fout_cov=open(sfcoverage,"wt")
        fout_cov.write(str(f))
        fout_cov.close()
    else:
        #read in file
        fin_cov=open(sfcoverage,"rt")
        f=float(fin_cov.readline())
        fin_cov.close()

    stest_output="Read coverage is: {0}".format(f)
    print stest_output
    return f

#calculate coverage
def calc_coverage(fpath, read_length, genome_length):
    cnt_lines=0
    cmd=""
    if fpath.lower().endswith(('.fq', '.fastq')):
        cmd="echo $(wc -l {0})".format(fpath)
    elif fpath.lower().endswith(('.fastq.gz', 'fq.gz', 'gz')):
        cmd="echo $(zcat {0} | wc -l)".format(fpath)
    else:
        print "Something wrong with the raw reads files format:", fpath

    print_command("Running command: "+ cmd +"...")

    tp=tuple(Popen(cmd, shell = True, stdout = PIPE).communicate())
    lcnt=str(tp[0]).split()
    print_command("The number of reads in file {0} is {1}".format(fpath,tp))
    cnt_lines=int(lcnt[0])
    cnt_reads=int(cnt_lines)/4
    cov=float(cnt_reads*read_length)/float(genome_length)
    return cov


'''
Description:
    Read in a fasta file
Parameters:
    sffa: path of fa file
    brstrip: whether strip '\n' or others can be striped by rstrip()
'''
def read_contig_fa(sffa,brstrip):
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

