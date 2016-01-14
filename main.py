__author__ = 'Chong Chu'
#!/usr/bin/env python

import sys
import os
import argparse
from subprocess import *

from Utility import *
from Assembly import *
from BasicInfoPaser import *
from ClassifyContigs import *
from MergeContigs import *
from FilterAndScaffold import *

############################################################################76 characters

############################################################################
###public values #######################
MIN_REPEAT_FREQ=1000
RELATIVE_FREQ_ON=1
RANGE_ASM_FREQ_DEC=3.0
RANGE_ASM_FREQ_GAP=0.6
K_MIN=10
K_MAX=60
K_INC=10
K_DFT=50
READ_LENGTH=135
GENOME_LENGTH=3209300000
READ_DEPTH=1
ASM_NODE_LENGTH_OFFSET=-1
MIN_CONTIG_LENGTH=100
COV_DIFF_CUTOFF=0.5
MIN_SUPPORT_PAIRS=20
MIN_FULLY_MAP_RATIO=0.2
TR_SIMILARITY=0.8 ##tandem repeats similarity between each two contigs, used as threshold to classify TR

IS_DUPLICATE_REPEATS=0.95
RANGE_ASM_FREQ_INC_TIMES=5 ##increase 5 times
RANGE_ASM_FREQ_DEC_TIMES=0.1 ##decrease 0.1 times
RM_DUP_BF_MERGE_CUTOFF=0.9
RM_DUP_AF_MERGE_CUTOFF=0.85

bpaired=True
sfleft_reads=""
sfright_reads=""
sfsingle_reads=""
file_list=[]

local_BWA_PATH="bwa"
local_SAMTOOLS_PATH="samtools"
local_REFINER_PATH="./TERefiner_1"
local_JELLYFISH_PATH=""
local_VELVET_PATH=""
local_CONTIGS_MERGER_PATH="./ContigsMerger"
local_THREADS=15
local_OUTPUT_FOLDER="./REPdenovo_Output/"
local_VERBOSE=0
#############################################################################

#####read in configuration file##############################################
def read_configfile(sfconfig):
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
        elif parts[0]=='TREADS':
            global local_TREADS
            local_TREADS=int(parts[1])
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
            if parts[1]!="GLOBAL":
                global local_BWA_PATH
                local_BWA_PATH=parts[1]
        elif parts[0]=="SAMTOOLS_PATH":
            if parts[1]!="GLOBAL":
                global local_SAMTOOLS_PATH
                local_SAMTOOLS_PATH=parts[1]
        elif parts[0]=="REFINER_PATH":
            if parts[1]!="GLOBAL":
                global local_REFINER_PATH
                local_REFINER_PATH=parts[1]
        elif parts[0]=='JELLYFISH_PATH':
            if parts[1]!="GLOBAL":
                global local_JELLYFISH_PATH
                local_JELLYFISH_PATH=parts[1]
        elif parts[0]=='VELVET_PATH':
            if parts[1]!="GLOBAL":
                global local_VELVET_PATH
                local_VELVET_PATH=parts[1]
        elif parts[0]=="CONTIGS_MERGER_PATH":
            if parts[1]!="GLOBAL":
                global local_CONTIGS_MERGER_PATH
                local_CONTIGS_MERGER_PATH=parts[1]
        elif parts[0]=="OUTPUT_FOLDER":
            global local_OUTPUT_FOLDER
            local_OUTPUT_FOLDER=parts[1]
        elif parts[0]=="READ_DEPTH":
            global READ_DEPTH
            READ_DEPTH=float(parts[1])
        elif parts[0]=="GENOME_LENGTH":
            global GENOME_LENGTH
            GENOME_LENGTH=int(parts[1])
        elif parts[0]=='VERBOSE':
            global local_VERBOSE
            local_VERBOSE=int(parts[1])
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

    ##check each path exist
    if local_BWA_PATH!="bwa" and os.path.exists(local_BWA_PATH)==False:
        info="BWA is not installed or provided path({0}) is wrong.".format(local_BWA_PATH)
        sys.exit(info)
    if local_SAMTOOLS_PATH!="samtools" and os.path.exists(local_SAMTOOLS_PATH)==False:
        info="SAMTOOLS is not installed or provided path({0}) is wrong.".format(local_SAMTOOLS_PATH)
        sys.exit(info)
    if local_JELLYFISH_PATH!="" and os.path.exists(local_JELLYFISH_PATH)==False:
        info="Jellyfish is not installed or provided path({0}) is wrong.".format(local_JELLYFISH_PATH)
        sys.exit(info)
    if local_VELVET_PATH!="" and os.path.exists(local_VELVET_PATH)==False:
        info="Velvet is not installed or provided path({0}) is wrong.".format(local_VELVET_PATH)
        sys.exit(info)
    if os.path.exists(local_REFINER_PATH)==False:
        info="{0} is not found !".format(local_REFINER_PATH)
        sys.exit(info)
    if os.path.exists(local_CONTIGS_MERGER_PATH)==False:
        info="{0} is not found !".format(local_CONTIGS_MERGER_PATH)
        sys.exit(info)

    set_parameters(local_BWA_PATH,local_SAMTOOLS_PATH,local_REFINER_PATH,\
                  local_JELLYFISH_PATH,local_VELVET_PATH,\
                  local_THREADS,local_OUTPUT_FOLDER,local_VERBOSE)



######read in raw reads files###########################################################################################
def read_rawreads_list(sfreads_list):
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


############main procedure of assembly and process contigs#####################
def rmTR_from_contigs(vtr):
    global  local_OUTPUT_FOLDER

    dtr={}
    for vrecords in vtr:
        for sname in vrecords:
            dtr[sname]=1
    sall=local_OUTPUT_FOLDER+"contigs.fa"
    shutil.copy2(sall, local_OUTPUT_FOLDER+"contigs_before_remove_TR.fa")
    dcontigs=read_contig_fa(sall,False)

    fnew_contigs=open(sall,"w")
    for key in dcontigs:
        if dtr.has_key(key):
            continue
        else:
            fnew_contigs.write(">"+key)
            fnew_contigs.write(dcontigs[key])
    fnew_contigs.close()


def clear_all():
    global  local_OUTPUT_FOLDER

    cmd="rm {0}reads_coverage.txt".format(local_OUTPUT_FOLDER)
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print_command(cmd)
    cmd="rm {0}dumped_*mers.txt".format(local_OUTPUT_FOLDER)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    print_command(cmd)
    print "Running command: "+ cmd +"..."
    Popen(cmd, shell = True, stdout = PIPE).communicate()


def usage():
    print 'Usage: python {0} -c Options -g configure_file -r raw_reads_list_file(or sam file)\n'.format(sys.argv[0]),
    print 'Options:\n',
    print '    Clean          Remove all the temporary files\n',
    print '    All            Run the whole pipeline\n',
    print '    Assembly       Only run the assembly part\n',
    print '    RmDup          Remove duplicate and contained ones\n',
    print '    Classify       Classfy tandem repeats from other types of repeats',
    print '    Scaffolding    Only run the scaffolding part, with the given contig file\n'
    print '    Analysis       Scaffolding and analysis with given bam files and contig file\n',
    print 'Example of running the whole pipeline: python main.py -c All -g config.txt -r ram_reads.txt'
########################################################################################################################


def get_args():
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
        description='Run the pipeline of REPdenovo')
    # Add arguments
    parser.add_argument(
        '-g', '--config', type=str, help='Configuration file name', required=True)
    parser.add_argument(
        '-r', '--reads', type=str, help='Raw reads file names', required=True)
    parser.add_argument(
        '-c', '--command', type=str, help='Specific command', required=True)

    args = parser.parse_args()
    sfconfig=args.config
    sfreads_list = args.reads
    scommand=args.command

    return sfconfig, sfreads_list, scommand

#####Main Procedure#####################################################################################################
def main_func(scommand,sfconfig,sfreads_list):
    global READ_DEPTH
    global MIN_REPEAT_FREQ
    global local_OUTPUT_FOLDER
    global bpaired,file_list, sfleft_reads, sfright_reads,sfsingle_reads
    global READ_LENGTH,GENOME_LENGTH,COV_DIFF_CUTOFF,MIN_SUPPORT_PAIRS
    global K_MIN, K_MAX, K_INC, ASM_NODE_LENGTH_OFFSET, TR_SIMILARITY
    global RANGE_ASM_FREQ_DEC_TIMES, RANGE_ASM_FREQ_INC_TIMES,MIN_CONTIG_LENGTH,RANGE_ASM_FREQ_DEC
    global local_THREADS, RM_DUP_BF_MERGE_CUTOFF, RM_DUP_AF_MERGE_CUTOFF

    assert os.path.exists(sfconfig),"configuration file is not found"
    read_configfile(sfconfig)

    if local_OUTPUT_FOLDER[-1]!='/':
        local_OUTPUT_FOLDER+="/"
    if not os.path.exists(local_OUTPUT_FOLDER):
        os.makedirs(local_OUTPUT_FOLDER)

    assert os.path.exists(sfreads_list),"raw reads list file is not found"
    read_rawreads_list(sfreads_list)

    if bpaired==True:
        sfleft_reads=file_list[0][0]
        sfright_reads=file_list[1][0]
        ilinsert_size=file_list[0][1] # left insert size
        irinsert_size=file_list[1][1] # right insert size
    else:
        sfsingle_reads=file_list[0][0]
        isinsert_size=file_list[0][1]

    if READ_DEPTH<=1:
        READ_DEPTH=calc_total_coverage(bpaired,sfleft_reads,sfright_reads,\
                                     sfsingle_reads,READ_LENGTH,GENOME_LENGTH)

    MIN_REPEAT_FREQ=int(MIN_REPEAT_FREQ*READ_DEPTH) ##use value that relative to coverage

    if scommand=="Scaffolding":
        align_read_to_contigs(file_list)
        scaffold(file_list,READ_LENGTH,COV_DIFF_CUTOFF, MIN_SUPPORT_PAIRS)
    elif scommand=="Analysis":
        scaffold_with_bam_list(file_list,READ_LENGTH,COV_DIFF_CUTOFF,MIN_SUPPORT_PAIRS)
    elif scommand=="Clean":
        clear_all()
    elif scommand=='Assembly':
        assembly(K_MIN, K_MAX, K_INC, MIN_REPEAT_FREQ, RANGE_ASM_FREQ_DEC_TIMES,\
                 RANGE_ASM_FREQ_INC_TIMES,READ_DEPTH,\
                 ASM_NODE_LENGTH_OFFSET, MIN_CONTIG_LENGTH,RANGE_ASM_FREQ_DEC,\
                 bpaired, sfleft_reads, sfright_reads, sfsingle_reads)
        merge_contigs(local_CONTIGS_MERGER_PATH,local_OUTPUT_FOLDER,local_THREADS, \
                     RM_DUP_BF_MERGE_CUTOFF, RM_DUP_AF_MERGE_CUTOFF)
    elif scommand=="Merging":
        merge_contigs(local_CONTIGS_MERGER_PATH,local_OUTPUT_FOLDER,local_THREADS, \
                     RM_DUP_BF_MERGE_CUTOFF, RM_DUP_AF_MERGE_CUTOFF)
    elif scommand=='Classify':
        vtr=classify_contigs(local_OUTPUT_FOLDER, K_MIN, K_MAX, K_INC, \
                            ASM_NODE_LENGTH_OFFSET, TR_SIMILARITY)
        #rmTR_from_contigs(vtr)
    elif scommand=="All":
        assembly(K_MIN, K_MAX, K_INC, MIN_REPEAT_FREQ, RANGE_ASM_FREQ_DEC_TIMES, \
                 RANGE_ASM_FREQ_INC_TIMES,READ_DEPTH,\
                 ASM_NODE_LENGTH_OFFSET, MIN_CONTIG_LENGTH,RANGE_ASM_FREQ_DEC,\
                 bpaired, sfleft_reads, sfright_reads, sfsingle_reads)
        merge_contigs(local_CONTIGS_MERGER_PATH,local_OUTPUT_FOLDER, \
                     local_THREADS, RM_DUP_BF_MERGE_CUTOFF, RM_DUP_AF_MERGE_CUTOFF)
        #vtr=classifyContigs(OUTPUT_FOLDER, K_MIN, K_MAX, K_INC, ASM_NODE_LENGTH_OFFSET, TR_SIMILARITY)
        #rmTR_from_contigs(vtr)
        align_read_to_contigs(file_list)
        scaffold(file_list,READ_LENGTH,COV_DIFF_CUTOFF, MIN_SUPPORT_PAIRS)
    elif scommand=="RmDup":
        rm_dup_contain(local_OUTPUT_FOLDER, RM_DUP_AF_MERGE_CUTOFF)
    else:
        print "Wrong parameters!!!\n"
        usage()


if __name__ == "__main__":
    if len(sys.argv) <= 3:
        usage()
        raise SystemExit

    #scommand=sys.argv[1]
    #sfconfig=sys.argv[2]
    #sfreads_list=sys.argv[3]
    sfconfig, sfreads_list, scommand = get_args()
    main_func(scommand,sfconfig,sfreads_list)
