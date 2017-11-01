#!/usr/bin/env python
import sys
import glob
import argparse 
from collections import defaultdict

import pysam
import numpy as np

__author__ = 'Thomas Hansen (tbh@mbg.au.dk)'
__version__ = '1.0.0'



bp = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'N' : 'N',
}

def get_complement (seq):

    return "".join([bp[x] for x in seq])[::-1]

    
class locus (object):

    def __init__(self, chr, start, end, name, strand, samples, read_count, add = 10):
       
        global agotron_id
        
        self.chr, self.start, self.end, self.name, self.strand = chr, start, end, name, strand
        self.samples = samples
        self.rc = read_count
                
        self.id = agotron_id
        self.p5end = []
        self.read_length = []
        self.p5density = defaultdict (int)
        self.read_density = defaultdict (int)
        self.read_density2 = defaultdict (lambda: defaultdict(int)) # used for coverage
        self.total_reads = 0
        self.bam_reads = defaultdict (int)
        self.mapq = defaultdict (list)
        
        self.add = add
        
    
    def add_read (self, read, bam, mapq_cutoff = 13):
        
        if (self.strand == "+" and read.is_reverse) or (self.strand == "-" and not read.is_reverse):
            return

        if read.mapq < args.mapq:
            return
        
        pos5p =  read.pos if self.strand == "+" else read.aend
                        
        self.p5end.append (pos5p)
        self.p5density[pos5p] += 1
        
        rlength = abs(read.aend - read.pos) - 1
        
        rpos = pos5p
        
        for i in range (read.pos, read.aend):
            self.read_density[i] += 1
            self.read_density2[rlength][i] += 1
        
        
        self.read_length.append (rlength)        
        self.total_reads += 1
        self.bam_reads[bam] += 1
    
    
    def get_seq (self):
    
        global fastafile
    
        if fastafile == None:
            return ""
        
         
        seq  = fastafile.fetch (self.chr, int(self.start) - self.add, int (self.end) + self.add)

        if self.strand == "-":
            seq = get_complement (seq)
        
        return seq
        
    def get_total_reads (self):
    
        return self.total_reads
    
        
    def get_percentiles (self):
                
        return [np.median (self.read_length), np.percentile (self.read_length, 5),  np.percentile (self.read_length, 95)] 
        
    def get_5p_heterogeneity (self): 
        
        max_reads, max_pos, total_reads = 0, -1, self.get_total_reads ()
                
        for pos in self.p5density:
            
            (max_reads, max_pos) = (self.p5density[pos], pos) if self.p5density[pos] > max_reads else (max_reads, max_pos)
                               
        return [float (max_reads) / float (total_reads), max_pos]

    def get_5p_position (self): #deprecated
      
        range_s = self.start   if self.strand == "+" else self.end-2
        range_e = self.start+2 if self.strand == "+" else self.end
        
        preads = 0
        
        for pos in range (range_s, range_e):
            
            preads += self.p5density[pos]
                               
        return [float (preads) / float (self.total_reads)]
 
       
    def get_inside (self):
    
        total_density = sum(self.read_density.values())
        overlapping_reads = 0
        
        for i in range (self.start, self.end):
            overlapping_reads += self.read_density[i]
        
        return float (overlapping_reads) / float (total_density)
        
        
    def enough_reads (self, th, ts):
    
        nts = int(ts) if ts.isdigit() else len(self.bam_reads)
        nth = 0
        for key in self.bam_reads:
            if float (self.bam_reads[key]) >= float (self.rc[key]) * float (th) / 1000000.0:
                nth += 1
                return True
    
        return nth >= nts
        
    
    def output (self):
   
        exp = [float (self.bam_reads[s]) * 1000000.0 / float (self.rc[s]) for s in self.samples]
        counts = [self.bam_reads[s] for s in self.samples]
    
        return [str (s) for s in ([self.name, self.get_inside()] + self.get_percentiles () + self.get_5p_heterogeneity ()) + [self.get_seq ()] + exp + counts]
        
    
    def output_coverage (self, file):

        if file == None:
            return
            
        for length in self.read_density2:
            
            for pos in self.read_density2[length]:
            
                print >>file, "\t".join (["locus_" + str (self.id), str(length), str(pos), str(self.read_density2[length][pos])])


parser = argparse.ArgumentParser(description='Intersect and analyse BED and BAM files.')

parser.add_argument("-g", dest="genome", type=str, help="Genome fasta, must be indexed", required=True)
parser.add_argument("-f", dest="files", nargs="+", type=str, help="Input bam-files. If empty all bam-files in current directoray are used")
parser.add_argument("-c", dest="coverage", type=str, default = "coverage.txt", help="Output coverage filename") 
parser.add_argument('-tr', dest='tr', type=float, default=5, help="Threshold for RPMM expression (default=5)") #RPMM: reads per mapped million
parser.add_argument('-ts', dest='ts', type=str, default="2", help="How many samples to meet RPMM expression threshold. For all samples; 'all'. (default = 2)")
parser.add_argument('-m', dest='max_outside', type=float, default=0.1, help="Tolerance for reads mapping partly outside locus")
parser.add_argument('-q', dest='mapq', type=int, default=13, help="MAPQ cutoff")
parser.add_argument('-a', dest='add_flank', type=int, default=10, help="Add <int> flanking nucleotides (up and downstream) to the loci sequence output (default=10)")
parser.add_argument('-v', '--version', action='version', version='%(prog)s '  + __version__)

args = parser.parse_args()

agotron_id = 1

### ----------------------------------------------------------
### Check argument
### ----------------------------------------------------------

samples = args.files

if not samples or len (samples) == 0:
    samples = glob.glob ("*.bam")

samples.sort()
    

read_count = {}
bam_entries = {}

fastafile = None

if args.genome != "":
    fastafile = pysam.Fastafile (args.genome)


### ----------------------------------------------------------
### First: Get total number of mapped reads in bam using pysam
### ----------------------------------------------------------
    
print >>sys.stderr, "Getting mapped readcount..."

bam_files = defaultdict (pysam.AlignmentFile)

for bam in samples:

    nmapped = 0
        
    for idx in pysam.idxstats(bam).split('\n'):
        
        idxc = idx.strip().split("\t")
        nmapped += int (idxc[2]) if len (idxc) > 2 else 0

    
    print >>sys.stderr, "{} mapped reads: {}".format (bam, nmapped)
        
    read_count[bam] = nmapped
    bam_files[bam] = pysam.AlignmentFile(bam, "rb")
    
    header = bam_files[bam].header
    
    for dict in header["SQ"]:
        
        bam_entries[dict["SN"]] = dict["LN"]


fcoverage = None
        
if args.coverage != "":
    fcoverage = open (args.coverage, "w")
    print >>fcoverage, "\t".join (["id", "length", "pos", "reads"])


### ------------------------------------
### print header
### ------------------------------------

header_bed = ["#chrom", "start", "end", "name", "score", "strand"]
header_analysis = ["host_name", "overlap", "median", "percentile_5th", "percentile_95th", "p5homogeneity", "p5position", "seq"]
header_exp = ["RPMM_" + s for s in samples] + ["reads_" + s for s in samples]
    
print "\t".join (header_bed + header_analysis + header_exp)



### ------------------------------------
### Loop through all entries in bed-file
### ------------------------------------

print >>sys.stderr, "Processing loci..."

dDone = {} # Only consider unique 5'ss

for line in sys.stdin: 

    cols = line.strip().split('\t')
    
    if len(cols) < 6:
        continue
    
    if not cols[1].isdigit () or not cols[2].isdigit (): # header
        continue
    
    chr, start, end, name, strand = cols[0], int(cols[1]), int (cols[2]), cols[3], cols[5]
    
    if chr not in bam_entries:
        continue
 
    pos = (chr, start if strand == "+" else end)
        
    if pos in dDone:
        continue
    
    dDone[pos] = 1
    
    ilocus = locus (chr, start, end, name, strand, samples, read_count)
       
    for bam in samples:
                
        for read in bam_files[bam].fetch (chr, start, end):    
            
            ilocus.add_read (read, bam, args.mapq)
           
        
    if not ilocus.enough_reads (args.tr, args.ts) or ilocus.get_inside () < 1 - args.max_outside: # Too few reads or only 
        continue
    
     
    print "\t".join (str (s) for s in [chr, start, end, "locus_" + str(agotron_id), ilocus.get_total_reads(), strand] + ilocus.output())
    
    ilocus.output_coverage (fcoverage)

    agotron_id += 1
    

if args.coverage != "":
    fcoverage.close()
