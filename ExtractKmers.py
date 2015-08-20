__author__ = 'Chong Chu'
import sys
import os

ffile=open(sys.argv[1])
cutoff=int(sys.argv[2])
seq=""
cnt=0
for line in ffile:
    if line[0]==">":
        print seq.rstrip(),cnt
        cnt=int(line[1:])
    else:
        seq=line

if seq!="" and cnt>=cutoff:
    print seq,cnt
ffile.close()