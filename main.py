
__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
import shutil
from subprocess import *
from KmerCount import *
from Assembly import *
from BasicInfoPaser import *
from FilterPEReads import *
from ClassifyContigs import *
from MergeContigs import *

##########################################################################################################
###public values #######################
MIN_REPEAT_FREQ=1000
RELATIVE_FREQ_ON=1
RANGE_ASM_FREQ_DEC=3.0
#RANGE_ASM_FREQ_INC=100
RANGE_ASM_FREQ_GAP=0.6
K_MIN=10
K_MAX=60
K_INC=10
K_DFT=50
READ_LENGTH=135
JELLYFISH_TREADS=5
GENOME_LENGTH=3209300000
ASM_NODE_LENGTH_OFFSET=-1
MIN_CONTIG_LENGTH=100
COV_DIFF_CUTOFF=0.5
MIN_SUPPORT_PAIRS=20
MIN_FULLY_MAP_RATIO=0.2
TR_SIMILARITY=0.8 ##tandem repeats similarity between each two contigs, used as threshold to classify TR
BWA_PATH="bwa"
SAMTOOLS_PATH="samtools"
REFINER_PATH="./TERefiner_1"
JELLYFISH_PATH=""
VELVET_PATH=""
BWA_THREADS=5
OUTPUT_FOLDER="./temp/"
READ_DEPTH=1
VERBOSE=0
IS_DUPLICATE_REPEATS=0.95
RANGE_ASM_FREQ_INC_TIMES=5 ##increase 5 times
RANGE_ASM_FREQ_DEC_TIMES=0.1 ##decrease 0.1 times
RM_DUP_BF_MERGE_CUTOFF=0.9
RM_DUP_AF_MERGE_CUTOFF=0.85

bpaired=True
sfleft_reads=""
sfright_reads=""
ilinsert_size=0 # left insert size
irinsert_size=0 # right insert size
ilsd_is=0 # left standard derivation of insert size
irsd_is=0 # right standard derivation of insert size
sfsingle_reads=""
isinsert_size=0 # single-end reads insert size
issd_is=0 #single-end reads standard derivation of insert size
file_list=[]
########################################################################################################

#####read in configuration file#########################################################################
def readConfigFile(sfconfig):
    fconfig=open(sfconfig,"r")
    sflines=fconfig.readlines()
    for line in sflines:
        parts=line.split()
        if len(parts) < 2:
            print "Wrong paramters in configuration file!!!"
            continue

        if parts[0]=='MIN_REPEAT_FREQ' :
           global MIN_REPEAT_FREQ
           MIN_REPEAT_FREQ=int(parts[1])
        elif parts[0]=='RANGE_ASM_FREQ_DEC':
            global RANGE_ASM_FREQ_DEC
            RANGE_ASM_FREQ_DEC=float(parts[1])
        elif parts[0]=='RANGE_ASM_FREQ_GAP':
            global RANGE_ASM_FREQ_GAP
            RANGE_ASM_FREQ_GAP=float(parts[1])
        elif parts[0]=='K_MIN':
            global K_MIN
            K_MIN=int(parts[1])
        elif parts[0]=='K_MAX':
            global K_MAX
            K_MAX=int(parts[1])
        elif parts[0]=='K_INC':
            global K_INC
            K_INC=int(parts[1])
        elif parts[0]=='K_DFT':
            global K_DFT
            K_DFT=int(parts[1])
        elif parts[0]=='READ_LENGTH':
            global READ_LENGTH
            READ_LENGTH=int(parts[1])
        elif parts[0]=='JELLYFISH_TREADS':
            global JELLYFISH_TREADS
            JELLYFISH_TREADS=int(parts[1])
        elif parts[0]=="ASM_NODE_LENGTH_OFFSET":
            global ASM_NODE_LENGTH_OFFSET
            ASM_NODE_LENGTH_OFFSET=int(parts[1])
        elif parts[0]=='COV_DIFF_CUTOFF':
            global COV_DIFF_CUTOFF
            COV_DIFF_CUTOFF=float(parts[1])
        elif parts[0]=='MIN_SUPPORT_PAIRS':
            global MIN_SUPPORT_PAIRS
            MIN_SUPPORT_PAIRS=int(parts[1])
        elif parts[0]=='MIN_FULLY_MAP_RATIO':
            global MIN_FULLY_MAP_RATIO
            MIN_FULLY_MAP_RATIO=float(parts[1])
        elif parts[0]=="TR_SIMILARITY":
            global TR_SIMILARITY
            TR_SIMILARITY=float(parts[1])
        elif parts[0]=='BWA_PATH':
            if parts[0]!="GLOBAL":
                global BWA_PATH
                BWA_PATH=parts[1]
        elif parts[0]=="SAMTOOLS_PATH":
            if parts[0]!="GLOBAL":
                global SAMTOOLS_PATH
                SAMTOOLS_PATH=parts[1]
        elif parts[0]=="REFINER_PATH":
            if parts[0]!="GLOBAL":
                global REFINER_PATH
                REFINER_PATH=parts[1]
        elif parts[0]=='JELLYFISH_PATH':
            if parts[1]!="GLOBAL":
                global JELLYFISH_PATH
                JELLYFISH_PATH=parts[1]
        elif parts[0]=='VELVET_PATH':
            if parts[1]!="GLOBAL":
                global VELVET_PATH
                VELVET_PATH=parts[1]
        elif parts[0]=="BWA_THREADS":
            global BWA_THREADS
            BWA_THREADS=int(parts[1])
        elif parts[0]=="OUTPUT_FOLDER":
            global OUTPUT_FOLDER
            OUTPUT_FOLDER=parts[1]
        elif parts[0]=="READ_DEPTH":
            global READ_DEPTH
            READ_DEPTH=float(parts[1])
        elif parts[0]=="GENOME_LENGTH":
            global GENOME_LENGTH
            GENOME_LENGTH=int(parts[1])
        elif parts[0]=='VERBOSE':
            global VERBOSE
            VERBOSE=int(parts[1])
        elif parts[0]=='IS_DUPLICATE_REPEATS':
            global IS_DUPLICATE_REPEATS
            IS_DUPLICATE_REPEATS=float(parts[1])
        elif parts[0]=="RANGE_ASM_FREQ_INC_TIMES":
            global RANGE_ASM_FREQ_INC_TIMES
            RANGE_ASM_FREQ_INC_TIMES=float(parts[1])
        elif parts[0]=="RANGE_ASM_FREQ_DEC_TIMES":
            global RANGE_ASM_FREQ_DEC_TIMES
            RANGE_ASM_FREQ_DEC_TIMES=float(parts[1])
        elif parts[0]=="RM_DUP_BF_MERGE_CUTOFF":
            global RM_DUP_BF_MERGE_CUTOFF
            RM_DUP_BF_MERGE_CUTOFF=float(parts[1])
        elif parts[0]=="RM_DUP_AF_MERGE_CUTOFF":
            global RM_DUP_AF_MERGE_CUTOFF
            RM_DUP_AF_MERGE_CUTOFF=float(parts[1])

    fconfig.close()

######read in raw reads files###################################################################################################################
def readRawReadsList(sfreads_list):
    freads=open(sfreads_list,'r')
    sflines=freads.readlines()
    i=0
    while i<len(sflines):
        parts=sflines[i].split()
        if parts[0]=="#":##ignore comments
            i += 1
            continue

        if len(parts)!=4:
            print "File(raw reads) list format is wrong!!!"
            exit()
        if int(parts[1])==-1: #single-end reads
            global bpaired
            bpaired=False
            sfsingle_reads=parts[0]
            isinsert_size=float(parts[2])
            issd_is=float(parts[3])
            file_list.append([sfsingle_reads,isinsert_size,issd_is])
        else: # paired-end reads
            #left
            sfleft_reads=parts[0]
            lgid=parts[1]
            ilinsert_size=float(parts[2])
            ilsd_is=float(parts[3])

            #read in right
            i+=1
            parts=sflines[i].split()
            if len(parts)!=4:
                print "File(Reads) list format is wrong!!!"
                exit()
            sfright_reads=parts[0]
            rgid=parts[1]
            irinsert_size=float(parts[2])
            irsd_is=float(parts[3])

            if lgid != rgid :
                print("File(Reads) list format is wrong, not a pair!!!")
                exit()
            file_list.append([sfleft_reads,ilinsert_size,ilsd_is])
            file_list.append([sfright_reads,irinsert_size,irsd_is])
        i += 1

    freads.close()

def printCommand(cmd):
    if VERBOSE!=0:
        print "Running command: "+cmd+" ..."
############main procedure of assembly and process contigs##############################################################

def calcTotalCoverage():
    print "Calculating average coverage...."
    sfcoverage=OUTPUT_FOLDER+"reads_coverage.txt"
    f=2.0
    if READ_DEPTH>1:
        f=READ_DEPTH
    else:
        if os.path.exists(sfcoverage)!=True:
            #print "Test: ", sfleft_reads###########################################################################################################3
            f=calcCoverage(sfleft_reads,READ_LENGTH, GENOME_LENGTH,VERBOSE) + calcCoverage(sfright_reads,READ_LENGTH, GENOME_LENGTH,VERBOSE)
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

def clearBeforeAssembly():
    cmd="rm -r {0}Asm_*".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print cmd
    cmd="rm {0}*.fa".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print cmd
    cmd="rm {0}*.fastq".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print cmd

def assembly():

    clearBeforeAssembly()

    temp_k=K_MIN
    last_temp_k=-1

    print "Calculating average coverage...."
    sfcoverage=OUTPUT_FOLDER+"reads_coverage.txt"
    f=2.0
    if READ_DEPTH>1:
        f=READ_DEPTH
    else:
        if os.path.exists(sfcoverage)!=True:
            f=calcCoverage(sfleft_reads,READ_LENGTH, GENOME_LENGTH,VERBOSE) + calcCoverage(sfright_reads,READ_LENGTH, GENOME_LENGTH,VERBOSE)
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
            cntKmer(JELLYFISH_PATH, temp_k, JELLYFISH_TREADS, sfleft_reads, sfright_reads, min_dump_cnt, kmer_path, kmer_jf,VERBOSE)
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


def rmTRFromContigs(vtr):
    dtr={}
    for vrecords in vtr:
        for sname in vrecords:
            dtr[sname]=1
    sall=OUTPUT_FOLDER+"contigs.fa"
    shutil.copy2(sall, OUTPUT_FOLDER+"contigs_before_remove_TR.fa")
    dcontigs=readContigFa(sall,False)

    fnew_contigs=open(sall,"w")
    for key in dcontigs:
        if dtr.has_key(key):
            continue
        else:
            fnew_contigs.write(">"+key)
            fnew_contigs.write(dcontigs[key])
    fnew_contigs.close()

def preprocessSam(sfsam):
    #filter out those unqualified pairs
    sfsam_temp="{0}_temp.sam".format(sfsam)
    mapq=30
    filterSam(sfsam, mapq,sfsam_temp)

    cmd="samtools view -h -S -b {0}_temp.sam > {1}.bam".format(sfsam,sfsam)
    #print "Running command: "+ cmd +"..."
    printCommand(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools sort {0}.bam {1}.sort".format(sfsam,sfsam)
    #print "Running command: "+ cmd +"..."
    printCommand(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="samtools index {0}.sort.bam".format(sfsam)
    #print "Running command: "+ cmd +"..."
    printCommand(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #remove useless files
    #os.remove(sfsam) ##temporarily leave it there
    os.remove(sfsam_temp)
    os.remove("{0}.bam".format(sfsam))

##################################################################################
#####Align reads back to contigs##################################################
def alignReadToContigs():
    sall=OUTPUT_FOLDER+"contigs.fa"
    cmd="samtools faidx {0}".format(sall)
    #print "Running command: "+ cmd +"..."
    printCommand(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #check whether already indexed, if not then create
    bwa_index_path="{0}/{1}.bwt".format(OUTPUT_FOLDER,sall)
    if os.path.exists(bwa_index_path)!=True:
        cmd="{0} index {1}".format(BWA_PATH,sall)
        #print "Running command: "+ cmd +"..."
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    nfiles=len(file_list)
    if nfiles % 2 !=0 :
        print "Something wrong with the file list, it's not paired!!!"
        return -1
    i=0
    j=0
    while i<nfiles:
        sfleft_reads=file_list[i][0]
        i+=1
        sfright_reads=file_list[i][0]
        i+=1
        sfsam="{0}_{1}.sam".format(sall,j)
        cmd="{0} mem -t {1} {2} {3} {4} > {5}".format(BWA_PATH,BWA_THREADS,sall,sfleft_reads,sfright_reads, sfsam)
        #print "Running command: "+ cmd +"..."
        if os.path.exists(sfsam)==False:
            printCommand(cmd)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        #filter out those unqualified pairs
        sfsam_temp="{0}_{1}_temp.sam".format(sall,j)
        mapq=30
        filterSam(sfsam, mapq,sfsam_temp)

        cmd="{0} view -h -S -b {1}_{2}_temp.sam > {3}_{4}.bam".format(SAMTOOLS_PATH,sall,j,sall,j)
        #print "Running command: "+ cmd +"..."
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="{0} sort {1}_{2}.bam {3}_{4}.sort".format(SAMTOOLS_PATH,sall,j,sall,j)
        #print "Running command: "+ cmd +"..."
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="{0} index {1}_{2}.sort.bam".format(SAMTOOLS_PATH,sall,j)
        #print "Running command: "+ cmd +"..."
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        j+=1

def scaffold():
    sall=OUTPUT_FOLDER+"contigs.fa"
    nfiles=len(file_list)
    if nfiles % 2 !=0 :
        print "Something wrong with the file list!!!"
        return -1
    i=0
    j=0

    while i<nfiles:
        ##################################################################################
        ######Remove some with low fully mapped reads
        #high_qulity_fa="high_quality_contigs.fa"
        #cmd="TERefiner_1 -C -b {0}.sort.bam -r {1} -c {2} -l {3} -o {4}".format(sall,sall,MIN_FULLY_MAP_RATIO,READ_LENGTH,high_qulity_fa)
        #Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##################################################################################
        #####Scaffolding##################################################################
        INSERT_SIZE=file_list[i][1]
        SD_INSERT_SIZE=file_list[j][2]

        fbam_cov="{0}contigs.fa_{1}.sam_for_coverage.sorted.bam".format(OUTPUT_FOLDER,j)
        cmd="{0} -L -b {1}_{2}.sort.bam -r {3} -l {4} -c {5} -t {6} -o {7}{8}_contig_pairs_info.txt -m {9} -d {10} -v {11}".format(REFINER_PATH,sall,j,sall, READ_LENGTH, COV_DIFF_CUTOFF, MIN_SUPPORT_PAIRS,OUTPUT_FOLDER,j,INSERT_SIZE,SD_INSERT_SIZE,fbam_cov)
        #print "Running command: "+ cmd +"..."
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##output concatenated ones
        cmd="{0} -S -r {1} -s {2}{3}_contig_pairs_info.txt -o {4}{5}_merged.fa".format(REFINER_PATH,sall,OUTPUT_FOLDER,j,OUTPUT_FOLDER,j)
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        i+=2
        j+=1

def scaffoldwithBamList():
    sall=OUTPUT_FOLDER+"contigs.fa"
    assert os.path.exists(sall)," file is not found"
    nfiles=len(file_list)
    j=0
    while j<nfiles:
        ##################################################################################
        ######Remove some with low fully mapped reads
        #high_qulity_fa="high_quality_contigs.fa"
        #cmd="TERefiner_1 -C -b {0}.sort.bam -r {1} -c {2} -l {3} -o {4}".format(sall,sall,MIN_FULLY_MAP_RATIO,READ_LENGTH,high_qulity_fa)
        #Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##################################################################################
        #####Scaffolding##################################################################
        #global INSERT_SIZE
        INSERT_SIZE=file_list[j][1]
        SD_INSERT_SIZE=file_list[j][2]
        sbam=file_list[j][0]

        if os.path.exists(sbam)==False:
            sbam=OUTPUT_FOLDER+sbam
        if os.path.exists(sbam)==False:
            print sbam," is not exist!!!"
            j+=1
            continue
        fbam_cov="{0}contigs.fa_{1}.sam_for_coverage.sorted.bam".format(OUTPUT_FOLDER,j)
        cmd="{0} -L -b {1} -r {2} -l {3} -c {4} -t {5} -o {6}{7}_contig_pairs_info.txt -m {8} -d {9} -v {10}".format(REFINER_PATH,sbam,sall, READ_LENGTH, COV_DIFF_CUTOFF, MIN_SUPPORT_PAIRS,OUTPUT_FOLDER,j,INSERT_SIZE,SD_INSERT_SIZE,fbam_cov)
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##output concatenated ones
        cmd="{0} -S -r {1} -s {3}{4}_contig_pairs_info.txt -o {5}{6}_merged.fa".format(REFINER_PATH,sall,OUTPUT_FOLDER,j,OUTPUT_FOLDER,j)
        printCommand(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        j+=1

def clearAll():
    cmd="rm {0}reads_coverage.txt".format(OUTPUT_FOLDER)
    #print "Running command: "+ cmd +"..."
    printCommand(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #print "Running command: "+ cmd +"..."
    printCommand(cmd)
    cmd="rm {0}dumped_*mers.txt".format(OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    #cmd="rm -r {0}Asm_*".format(OUTPUT_FOLDER)
    printCommand(cmd)
    print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()

def usage():
    print 'Usage: python {0} Options configure_file raw_reads_list_file(or sam file)\n'.format(sys.argv[0]),
    print 'Options:\n',
    print '    Clean          Remove all the temporary files\n',
    print '    All            Run the whole pipeline\n',
    print '    Assembly       Only run the assembly part\n',
    print '    Classify       Classfy tandem repeats from other types of repeats',
    print '    Scaffolding    Only run the scaffolding part, with the given contig file\n',
    print '    Preprocess     Preprocess sam file, and generate sorted and indexed bam file\n',
    print '    Analysis       Scaffolding and analysis with given bam files and contig file\n',
    print 'Example of running the whole pipeline: python main.py All configre.txt ram_reads.txt'
########################################################################################################################
#####Main Procedure#####################################################################################################

if len(sys.argv) <= 3:
    usage()
    raise SystemExit

sfconfig=sys.argv[2]
assert os.path.exists(sfconfig),"configuration file is not found"
readConfigFile(sfconfig)

if OUTPUT_FOLDER[-1]!='/':
    OUTPUT_FOLDER+="/"
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

if sys.argv[1]=="Preprocess":
    preprocessSam(sys.argv[3])
else:
    sfreads_list=sys.argv[3]
    assert os.path.exists(sfreads_list),"raw reads list file is not found"
    readRawReadsList(sfreads_list)

    if bpaired==True:
        sfleft_reads=file_list[0][0]
        sfright_reads=file_list[1][0]
        ilinsert_size=file_list[0][1] # left insert size
        irinsert_size=file_list[1][1] # right insert size
    else:
        sfsingle_reads=file_list[0][0]
        isinsert_size=file_list[0][1]

    if READ_DEPTH==1:
        READ_DEPTH=calcTotalCoverage()

    MIN_REPEAT_FREQ=int(MIN_REPEAT_FREQ*READ_DEPTH) ##use value that relative to coverage

    if sys.argv[1]=="Scaffolding":
        alignReadToContigs()
        scaffold()
    elif sys.argv[1]=="Analysis":
        scaffoldwithBamList()
    elif sys.argv[1]=="Clean":
        clearAll()
    elif sys.argv[1]=='Assembly':
        assembly()
        mergeContigs(OUTPUT_FOLDER, RM_DUP_BF_MERGE_CUTOFF, RM_DUP_AF_MERGE_CUTOFF)

    elif sys.argv[1]=='Classify':
        vtr=classifyContigs(OUTPUT_FOLDER, K_MIN, K_MAX, K_INC, ASM_NODE_LENGTH_OFFSET, TR_SIMILARITY)
        #rmTRFromContigs(vtr)
    elif sys.argv[1]=="All":
        assembly()
        mergeContigs(OUTPUT_FOLDER, RM_DUP_BF_MERGE_CUTOFF, RM_DUP_AF_MERGE_CUTOFF)
        #vtr=classifyContigs(OUTPUT_FOLDER, K_MIN, K_MAX, K_INC, ASM_NODE_LENGTH_OFFSET, TR_SIMILARITY)
        #rmTRFromContigs(vtr)
        alignReadToContigs()
        scaffold()
    else:
        print "Wrong parameters!!!\n"
        usage()

