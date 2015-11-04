__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
import shutil
from subprocess import *
from KmerCount import *
from Utility import OUTPUT_FOLDER
from Utility import getOutputFolder
from Utility import VELVET_PATH
from Utility import getVelvetPath
from Utility import JELLYFISH_PATH
from Utility import getJellyfishPath
from Utility import VERBOSE
from Utility import getVerbose
from Utility import THREADS
from Utility import getThreadsNum

def assembly(K_MIN, K_MAX, K_INC, MIN_REPEAT_FREQ, RANGE_ASM_FREQ_DEC_TIMES, RANGE_ASM_FREQ_INC_TIMES,READ_DEPTH,\
             ASM_NODE_LENGTH_OFFSET, MIN_CONTIG_LENGTH,RANGE_ASM_FREQ_DEC,\
             bpaired, sfleft_reads, sfright_reads, sfsingle_reads):

    OUTPUT_FOLDER=getOutputFolder()
    THREADS=getThreadsNum()
    VERBOSE=getVerbose()
    JELLYFISH_PATH=getJellyfishPath()
    VELVET_PATH=getVelvetPath()

    clearBeforeAssembly()
    temp_k=K_MIN
    last_temp_k=-1

    f=READ_DEPTH
    min_dump_cnt=0
    if int(f)>MIN_REPEAT_FREQ:
        min_dump_cnt=int(f)
    else:
        min_dump_cnt=int(MIN_REPEAT_FREQ)

    while temp_k <= K_MAX:
        sprnt="Working on {0}mer now...".format(temp_k)
        print sprnt

        ##kmer counting
        kmer_path='dumped_{0}mers.txt'.format(temp_k)
        kmer_path=OUTPUT_FOLDER+kmer_path
        kmer_jf=OUTPUT_FOLDER+"mer_counts.jf"
        if os.path.exists(kmer_path)!=True:
            if bpaired==True:
                cntKmer(JELLYFISH_PATH, temp_k, THREADS, sfleft_reads, sfright_reads, min_dump_cnt, kmer_path, kmer_jf,VERBOSE)
            else:
                cntKmer(JELLYFISH_PATH, temp_k, THREADS, sfsingle_reads, "-1", min_dump_cnt, kmer_path, kmer_jf,VERBOSE)

            if os.path.exists(kmer_jf) == True:
                os.remove(kmer_jf)

        ##Load all the kmers into memory
        dallkmers={}
        akmer_freq=0
        max_freq=0
        fdumped_kmer=open(kmer_path,'r')
        for akmer in fdumped_kmer:
            if akmer[0]=='>':
                akmer_freq=int(akmer[1:])
                if akmer_freq>max_freq:
                    max_freq=akmer_freq
            else:
                if akmer_freq>=int(MIN_REPEAT_FREQ):
                    dallkmers[akmer]=akmer_freq
        fdumped_kmer.close()

        print "Highest {0}mer frequency is {1}".format(temp_k,max_freq)

        cur_kmer_freq=max_freq
        cur_f1=0
        cur_f2=0
        output_folder1=''
        scontig_k="contigs_{0}mer.fa".format(temp_k)
        scontig_k=OUTPUT_FOLDER+scontig_k
        sconcate="{0}mer.temp_contigs.fa".format(temp_k)
        sconcate=OUTPUT_FOLDER+sconcate

        iround=1
        while cur_kmer_freq>=MIN_REPEAT_FREQ:
            cur_f1=int(cur_kmer_freq*RANGE_ASM_FREQ_DEC_TIMES)
            cur_f2=int(cur_kmer_freq*RANGE_ASM_FREQ_INC_TIMES)

            if cur_f1<MIN_REPEAT_FREQ:
                cur_f1=MIN_REPEAT_FREQ
            if cur_f2>max_freq:
                cur_f2=max_freq

            sss="{0}mer frequency are:{1} and {2}".format(temp_k,cur_f1,cur_f2)
            print sss

            kmer_fastq=OUTPUT_FOLDER+"kmers_fq.fastq"

            fout_fq=open(kmer_fastq,"w")
            ifq_cnt=0

            for kmer_seq in dallkmers:
                if dallkmers[kmer_seq]>cur_f2:
                    #dallkmers.pop(kmer_seq,None)
                    donothing=0
                elif dallkmers[kmer_seq]>=cur_f1 and dallkmers[kmer_seq]<=cur_f2:
                    #write into file
                    fout_fq.write("@seq{0}\n".format(ifq_cnt))
                    ifq_cnt=ifq_cnt+1
                    fout_fq.write(kmer_seq+"\n")
                    fout_fq.write("+\n")
                    seq_len=len(kmer_seq)
                    iss=1
                    quality=""
                    while iss<seq_len:
                        quality+="5"
                        iss+=1
                    fout_fq.write(quality+"\n")

            fout_fq.close()

            velvet_k=temp_k+ASM_NODE_LENGTH_OFFSET

            output_folder1="{0}Asm_{1}_{2}_{3}".format(OUTPUT_FOLDER,temp_k,cur_kmer_freq,velvet_k)
            sexist=output_folder1+"/contigs.fa"
            if os.path.exists(sexist)!=True:
                assembleKmer(VELVET_PATH,velvet_k, kmer_fastq, MIN_CONTIG_LENGTH, output_folder1,VERBOSE)

            outfile=open(sconcate, 'a')
            readfile1=open(sexist, 'r')
            #shutil.copyfileobj(readfile1, outfile)
            for aline in readfile1:
                if aline[0]==">":
                    aline=aline.strip()+"_{0}_{1}_{2}\n".format(temp_k,int(cur_f1),int(cur_f2))
                outfile.write(aline)
            readfile1.close()
            outfile.close()

            cur_kmer_freq=int(float(cur_kmer_freq)/RANGE_ASM_FREQ_DEC)
            iround=iround+1

        sall=OUTPUT_FOLDER+"contigs.fa"
        if K_MIN==K_MAX:
            scontig_backup=OUTPUT_FOLDER+"contigs.fa"
            shutil.copy2(sconcate,scontig_backup)
            print "kmin and kmax are equal....."

        with open(sall,"a") as fout:
            with open(sconcate,"r") as fin:
                for line in fin:
                    fout.write(line)

        last_temp_k=temp_k
        temp_k += K_INC

def clearBeforeAssembly():
    OUTPUT_FOLDER=getOutputFolder()
    cmd="rm -r {0}Asm_*".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print cmd
    cmd="rm {0}*.fa".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print cmd
    cmd="rm {0}*.fastq".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print cmd

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


def assembleKmer(vpath, k_len, ffastq, min_contig_len, output_folder, verbose):
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
    if verbose != 0:
        print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} {1} -min_contig_lgth {2} -unused_reads yes".format(vpath_g,newpath,min_contig_len)
    if verbose != 0:
        print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()
