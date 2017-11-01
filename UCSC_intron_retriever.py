#!/usr/bin/env python

import sys
import argparse

import MySQLdb as mdb


class sql (object):

    def __init__(self, cur):
        
        self.c2i = self.col2idx (cur)        
    
    def get (self, row, column):

        return row[self.c2i[column]]
        
    def geti (self, row, column):
    
        return (int(self.get(row, column)))
    
    def col2idx (self, cur):

        columns = cur.description 
       
        result = [column[0] for column in columns] 
        index = 0
        col = {}
        for r in result:
           col[r] = index
           index+=1

        return col
    
    
   
parser = argparse.ArgumentParser(description='Get intron coordinates from UCSC mySQL database.')

parser.add_argument('-db', dest='database', default="hg19", type=str, help="MySQL database (default='hg19')")
parser.add_argument('-table', dest='table', default="refGene", type=str, help="MySQL table (default='refGene')")
parser.add_argument('-max', dest='max_intron_length', default=150, type=int, help="Maximum intron length (default=150)")
parser.add_argument('-min', dest='min_intron_length', default=50, type=int, help="Minimum intron length (default=50)")

args = parser.parse_args()

db = args.database
table = args.table


con = None

dDone = {}

try:

    # connect to UCSC database
    con = mdb.connect('genome-mysql.cse.ucsc.edu', 'genome', '', db)
   
    with con:
           
        cur = con.cursor()
        
        query = "SELECT * FROM {} WHERE exonCount > 1".format (table)
               
        cur.execute(query)
        
                
        isql = sql(cur)
                
        row = cur.fetchall()
                
        for r in row:
            
            chr, start, end, strand, name = isql.get (r,"chrom"), isql.geti (r,"txStart"), isql.geti (r,"txEnd"), isql.get (r, "strand"), isql.get (r, "name2")             
            es = isql.get (r, "exonStarts").strip().split (",")
            ee = isql.get (r, "exonEnds").strip().split (",")
            ne = isql.get (r, "exonCount")
                          
            # If agotrons, extract short introns
           
            
            for i in range (1, ne):
                                
                intronStart = int(ee[i-1])
                intronEnd = int(es[i])
                
                iintron = i if strand == "+" else ne-i
                
                if abs (intronEnd - intronStart) <= args.max_intron_length and \
                   abs (intronEnd - intronStart) >= args.min_intron_length and \
                   not (chr, intronStart, intronEnd, strand) in dDone:
                                            
                    bedname = name
                    
                    print "\t".join ([chr, str (intronStart), str (intronEnd), bedname, "-1", strand])
                    
                    dDone[(chr, intronStart, intronEnd, strand)] = 1
         
            
            
        print >>sys.stderr, "Found", len(dDone), "unique coordinates"
    
except mdb.Error, e:
  
    print "Error %d: %s" % (e.args[0], e.args[1])
    sys.exit(1)

finally:
    
    if con:
        con.close()

        
        
