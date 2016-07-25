__author__ = 'Chong Chu'

import sys
import os
from subprocess import *
from Utility import print_command
from Utility import OUTPUT_FOLDER
from Utility import get_output_folder
from Utility import BWA_PATH
from Utility import get_bwa_path
from Utility import THREADS
from Utility import get_threads_num
from Utility import SAMTOOLS_PATH
from Utility import get_samtools_path
from Utility import REFINER_PATH
from Utility import get_refiner_path
from FilterPEReads import filter_sam

def preprocess_sam(sfsam):
    SAMTOOLS_PATH=get_samtools_path()
    #filter out those unqualified pairs
    sfsam_temp="{0}_temp.sam".format(sfsam)
    mapq=30
    filter_sam(sfsam, mapq,sfsam_temp)

    cmd="{0} view -h -S -b {1}_temp.sam > {2}.bam".format(SAMTOOLS_PATH,sfsam,sfsam)
    #print "Running command: "+ cmd +"..."
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} sort {1}.bam -o {2}.sort.bam".format(SAMTOOLS_PATH,sfsam,sfsam)
    #print "Running command: "+ cmd +"..."
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()
    cmd="{0} index {1}.sort.bam".format(SAMTOOLS_PATH,sfsam)
    #print "Running command: "+ cmd +"..."
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #remove useless files
    #os.remove(sfsam) ##temporarily leave it there
    os.remove(sfsam_temp)
    os.remove("{0}.bam".format(sfsam))

##################################################################################
#####Align reads back to contigs##################################################
def align_read_to_contigs(file_list):
    OUTPUT_FOLDER=get_output_folder()
    BWA_PATH=get_bwa_path()
    THREADS=get_threads_num()
    SAMTOOLS_PATH=get_samtools_path()

    sall=OUTPUT_FOLDER+"contigs.fa"
    cmd="{0} faidx {1}".format(SAMTOOLS_PATH,sall)
    #print "Running command: "+ cmd +"..."
    print_command(cmd)
    Popen(cmd, shell = True, stdout = PIPE).communicate()

    #check whether already indexed, if not then create
    bwa_index_path="{0}/{1}.bwt".format(OUTPUT_FOLDER,sall)
    if os.path.exists(bwa_index_path)!=True:
        cmd="{0} index {1}".format(BWA_PATH,sall)
        #print "Running command: "+ cmd +"..."
        print_command(cmd)
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
        cmd="{0} mem -t {1} {2} {3} {4} > {5}".format(BWA_PATH,THREADS,sall,sfleft_reads,sfright_reads, sfsam)
        #print "Running command: "+ cmd +"..."
        if os.path.exists(sfsam)==False:
            print_command(cmd)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        #filter out those unqualified pairs
        sfsam_temp="{0}_{1}_temp.sam".format(sall,j)
        mapq=30
        filter_sam(sfsam, mapq,sfsam_temp)

        cmd="{0} view -h -S -b {1}_{2}_temp.sam > {3}_{4}.bam".format(SAMTOOLS_PATH,sall,j,sall,j)
        #print "Running command: "+ cmd +"..."
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="{0} sort {1}_{2}.bam -o {3}_{4}.sort.bam".format(SAMTOOLS_PATH,sall,j,sall,j)
        #print "Running command: "+ cmd +"..."
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        cmd="{0} index {1}_{2}.sort.bam".format(SAMTOOLS_PATH,sall,j)
        #print "Running command: "+ cmd +"..."
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        j+=1

def scaffold(file_list,READ_LENGTH,COV_DIFF_CUTOFF, MIN_SUPPORT_PAIRS):
    OUTPUT_FOLDER=get_output_folder()
    REFINER_PATH=get_refiner_path()

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
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##output concatenated ones
        cmd="{0} -S -r {1} -s {2}{3}_contig_pairs_info.txt -o {4}{5}_merged.fa".format(REFINER_PATH,sall,OUTPUT_FOLDER,j,OUTPUT_FOLDER,j)
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        i+=2
        j+=1

def scaffold_with_bam_list(file_list,READ_LENGTH,COV_DIFF_CUTOFF,MIN_SUPPORT_PAIRS):
    OUTPUT_FOLDER=get_output_folder()
    REFINER_PATH=get_refiner_path()

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
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        ##output concatenated ones
        cmd="{0} -S -r {1} -s {3}{4}_contig_pairs_info.txt -o {5}{6}_merged.fa".format(REFINER_PATH,sall,OUTPUT_FOLDER,j,OUTPUT_FOLDER,j)
        print_command(cmd)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

        j+=1
