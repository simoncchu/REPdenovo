from MergeContigs import *

fcontig=sys.argv[1]
foutput=sys.argv[2]
cutoff=0.95
brm_contained=True
removeDuplicateContained(fcontig, foutput, cutoff, brm_contained)