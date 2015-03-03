__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
from subprocess import *

#jellyfish count -m 35 -s 500M --bf-size 30M -t 5 -C -F 2 <(zcat NA12878_ERR194147_1.fastq.gz) <(zcat NA12878_ERR194147_2.fastq.gz)
def cntKmer(jpath, k_len, ithreads, flreads, frreads, min_cnt, foutput_dump, foutput_jf, VERBOSE):
    print('Counting kmers...')

    if jpath!="" and jpath[-1]!="/":
        jpath+="/"
    jpath+="jellyfish"

    if jpath != "jellyfish":
        assert os.path.exists(jpath),"jellyfish is not installed, or given the wrong path in configuration file!!!"

    #check file format is fq or fastq.gz
    cmd=""
    if flreads.lower().endswith(('.fq', '.fastq')) and frreads.lower().endswith(('.fq', '.fastq')):
        cmd='{0} count -s 500M --bf-size 30M -C -m {1} -t {2} -o {3} -F 2 {4} {5}'.format(jpath, k_len, ithreads, foutput_jf, flreads, frreads)
        if VERBOSE != 0:
            print "Running command: "+ cmd +"..."
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    elif flreads.lower().endswith(('.fastq.gz', 'fq.gz', 'gz')):
        cmd_temp='ls {0} {1} | xargs -n 1 echo gunzip -c > {2}_generators'.format(flreads, frreads, foutput_dump)
        if VERBOSE != 0:
            print "Running command: "+ cmd_temp +"..."
        Popen(cmd_temp, shell = True, stdout = PIPE).communicate()
        cmd = '{0} count -s 500M --bf-size 30M -C -m {1} -t {2} -o {3} -F 2 -g {4}_generators -G 2'.format(jpath, k_len, ithreads, foutput_jf, foutput_dump)
        if VERBOSE != 0:
            print "Running command: "+ cmd +"..."
        Popen(cmd, shell = True, stdout = PIPE).communicate()
    else:
        print "The reads file is not ended with .fq, .fastq, .fq.gz, .fastq.gz or .gz!!!!"


    #dump
    cmd="{0} dump -L 50 -o {2} {3}".format(jpath, min_cnt,foutput_dump, foutput_jf)
    if VERBOSE != 0:
        print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()



