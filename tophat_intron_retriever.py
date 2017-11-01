#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict
    

### -----------------------------------------------------------------------------------------------------
### Get unique intron coordinates from tophat output, junctions.bed, to use for agotron detection 
###
### USAGE: python tophat_intron_retriever.py [ARGS] < junctions.bed > introns.bed
### -----------------------------------------------------------------------------------------------------
    
   
parser = argparse.ArgumentParser(description='Get intron coordinates from tophat output: junctions.bed')
parser.add_argument('-max', dest='max_intron_length', default=150, type=int)
parser.add_argument('-min', dest='min_intron_length', default=50, type=int)

args = parser.parse_args()

dReads = defaultdict (int)

for line in sys.stdin:

    cols = line.strip().split("\t")

    if len(cols) < 10:
        continue
        
    if not cols[1].isdigit () or not cols[2].isdigit ():
        continue
    
    chr, start, end, name, score, strand, junc = cols[0], int(cols[1]), int(cols[2]), cols[3], int(cols[4]), cols[5], cols[10].split(",")
    
    junc_s = int(junc[0])
    junc_e = int(junc[1])
    
    intron_start = start+junc_s
    intron_end = end-junc_e
    
    intron_length = abs (intron_end - intron_start)
    
    if intron_length > args.min_intron_length and intron_length < args.max_intron_length:        
        dReads[(chr, start+junc_s, end-junc_e, strand)] += score 
    

iintron = 1
    
for (chr, s, e, strand) in dReads:

    print "\t".join ([chr, str(s), str(e), "intron" + str(iintron), str (dReads[(chr, s, e, strand)]), strand])
    
    iintron += 1
    
    
print >>sys.stderr, "Found", len (dReads), "unique coordinates"