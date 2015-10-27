__author__ = 'Chong Chu'
import sys
import os
import shutil
from subprocess import *
from BasicInfoPaser import *

if len(sys.argv) < 2:
    print "Please check parameter usage!!!"
    raise SystemExit

sfcontig=sys.argv[1]
assert os.path.exists(sfcontig),"contig file is not found"

###remove directly duplicates
sfcontig_no_dir_dup=sfcontig+"_no_direct_duplicate.fa"
READ_LENGTH=100

cmd="sh removeRepeatsOfOneContigSet.sh {0} {1} {2}".format(sfcontig, sfcontig_no_dir_dup, READ_LENGTH)
Popen(cmd, shell = True, stdout = PIPE).communicate()

SAME_CUTOFF=float(sys.argv[2])

###remove rotate duplicates 
sfcontig_rotate_same=sfcontig+"_rotate_supplementary_reverse_similar_contigs.txt"
frotate=open(sfcontig_rotate_same,"w")
dcontig=readContigFa(sfcontig_no_dir_dup,True)
for key1 in dcontig:
    seq1=dcontig[key1]
    for key2 in dcontig:
        if key1==key2: continue
        seq2=dcontig[key2]

        cmd="./TERefiner_1 -A -r {0} -s {1}".format(seq1,seq2)
        tp=tuple(Popen(cmd, shell = True, stdout = PIPE).communicate())
        parts=tp[0].split()
        ilen=int(parts[0])
        min_len=0
        if len(seq1)>len(seq2):
            min_len=len(seq2)
        else:
            min_len=len(seq1)

        if ilen > int(min_len*SAME_CUTOFF):
            frotate.write(key1+" ")
            frotate.write(key2+" \n")
frotate.close()