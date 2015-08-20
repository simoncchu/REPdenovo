__author__ = 'Chong Chu'

import sys
import os
from subprocess import *
from collections import defaultdict
from BasicInfoPaser import readContigFa


def isSufPreOverlap(seq,ilen,icutoff):
    if len(seq)<ilen:
        ilen=len(seq)

    sprefix=seq[0:ilen]
    istart=-1*ilen
    ssuffix=seq[istart:]

    cmd="./TERefiner_1 -M -r {0} -s {1}".format(sprefix,ssuffix)
    tp=tuple(Popen(cmd, shell = True, stdout = PIPE).communicate())
    parts=tp[0].split()
    pstart=int(parts[0])
    pend=int(parts[1])
    sstart=int(parts[2])
    send=int(parts[3])
    ioverlap=pend-pstart+1

    if ioverlap>=icutoff and pstart==1 and send==ilen:
        return True
    else:
        return False

'''
#parse out those suffix prefix overlapped contigs
def parseSPRepeatContig(sffa):
    dcontigs={}
    name=""
    seq=""
    last_name=""
    ffa=open(sffa)
    for line in ffa:
        line=line.rstrip()
        if line[0]=='>':
            name=line[1:]
            if last_name!="":
                if isSufPreOverlap(seq,40,7)==True:
                    dcontigs[last_name]=seq
                    #print last_name, seq ####################################################################
                seq=""
            last_name=name
        else:
            seq=seq+line
    ##check whether there is suffix and prefix overlap
    if isSufPreOverlap(seq,40,7)==True:
        dcontigs[last_name]=seq
        print last_name, seq ################################################################################
    ffa.close()
    return dcontigs
'''

def classifyContigs(OUTPUT_FOLDER, k_start, k_end, k_inc, ASM_NODE_LENGTH_OFFSET, TR_SIMILARITY):
    k_rounds=(k_end-k_start)/k_inc + 1
    #first read in all the length in to an array
    dseqs=[] #[]{}
    contigs_lenth=[] ##[list]{dict}{list}
    k_begin=k_start
    while k_begin<=k_end:
        sfcontig=OUTPUT_FOLDER+"contig_backup_for_{0}mer.fa".format(k_begin)

        dfa=readContigFa(sfcontig,True)
        dseqs.append(dfa)

        cmd="samtools faidx {0}".format(sfcontig)
        Popen(cmd, shell = True, stdout = PIPE).communicate()
        sfai=sfcontig+".fai"
        assert os.path.exists(sfai)," file is not found"
        sffai=open(sfai,"r")

        cntg_lenth=defaultdict(list)
        for line in sffai:
            parts=line.split()
            cntg_lenth[int(parts[1])].append(parts[0])

        contigs_lenth.append(cntg_lenth)
        k_begin=k_begin+k_inc
        sffai.close()

    #check the increase of length
    if len(contigs_lenth)<=0:
        print "No data for classification!!!!"
        return

    #find out all candidate groups
    vcandidates=[] ##
    sfrepeats=OUTPUT_FOLDER+"candidate_repeats.txt"
    frepeats=open(sfrepeats,"w")
    for key in contigs_lenth[0]:
        bRpt=True
        temp_key=key
        for i in range(1,k_rounds):
            temp_key=temp_key+k_inc
            if(contigs_lenth[i].has_key(temp_key)==False):
                bRpt=False
                break

        vrecord=[]
        if bRpt==True:
            stemp=""
            temp_key=key
            for m in range(0,k_rounds):
                vround=[]
                for scntg in contigs_lenth[m][temp_key]:
                    stemp=stemp+scntg+","
                    vround.append(scntg)
                vrecord.append(vround)
                stemp=stemp+" "
                temp_key=temp_key+k_inc
            stemp=stemp+"\n"+"\n"
            frepeats.write(stemp)

        if len(vrecord)>0:
            vcandidates.append(vrecord)
    frepeats.close()

    if len(vcandidates)<=0:
        return

    #print vcandidates #################################################################################################

    sftr=OUTPUT_FOLDER+"tandem_repeats_list.txt"
    ftr=open(sftr,"w")
    vtr=[]
    for i in range(0,len(vcandidates)): ###iterate each record
        '''
        #concatenate round[1:] to seperate strings
        vscon=[]
        ssep="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
        for j in range(1,k_rounds):### for each round
            sscon_temp=""
            for node_namei in vcandidates[i][j]:
                sscon_temp=sscon_temp+dseqs[j][node_namei]
                sscon_temp=sscon_temp+ssep
            vscon.append(sscon_temp)
        '''
        #for round 0, check each seq whetcher find a similar alignment with strings in vscon[]
        for node_name0 in vcandidates[i][0]:
            seq0=dseqs[0][node_name0]
            btr=True ##flag, is tandem repeat or not

            imin_len=k_start+ASM_NODE_LENGTH_OFFSET-1
            if isSufPreOverlap(seq0,k_end+1,imin_len)==False: #######################################suffix prefix overlap
                btr=False
                continue
            ##check whether can find an alignment in every other rounds
            voptimal_cntgs=[]

            for j in range(1,k_rounds):### for each round alignment
                max_len=0
                opt_node=""
                for node_namei in vcandidates[i][j]:
                    cmd="./TERefiner_1 -A -r {0} -s {1}".format(seq0,dseqs[j][node_namei])
                    tp=tuple(Popen(cmd, shell = True, stdout = PIPE).communicate())
                    parts=tp[0].split()
                    #istart=int(parts[0])
                    #iend=int(parts[1])
                    #ilen=iend-istart
                    ilen=int(parts[0])
                    if ilen>max_len:
                        max_len=ilen
                        opt_node=node_namei

                #print "max_len ",max_len#############################################################################################################
                #print "threshold ",int((len(seq0)-k_start)*TR_SIMILARITY)############################################################################

                if max_len > int((len(seq0)-k_start)*TR_SIMILARITY):
                    voptimal_cntgs.append(opt_node)
                else:
                    btr=False
                    break

                imin_len=j*k_inc+k_start+ASM_NODE_LENGTH_OFFSET-1
                if isSufPreOverlap(dseqs[j][opt_node],k_end+1,imin_len)==False: #######################################suffix prefix overlap
                    btr=False
                    break

            if btr==True:
                vtr.append(voptimal_cntgs)
                sout=">"+node_name0+"&"
                for opt_cntg in voptimal_cntgs:
                    sout=sout+opt_cntg+"_"
                sout=sout+str(len(voptimal_cntgs))
                sout=sout+"\n"
                ftr.write(sout)
                itr_len=len(seq0)-(k_start+ASM_NODE_LENGTH_OFFSET-1)
                tr=seq0[0:itr_len]
                ipos_seq=0
                tr_new=""
                while(ipos_seq<itr_len):
                    tr_new=tr_new+tr[ipos_seq:ipos_seq+60]
                    tr_new=tr_new+"\n"
                    ipos_seq+=60
                ftr.write(tr_new)

    ftr.close()
    return vtr
