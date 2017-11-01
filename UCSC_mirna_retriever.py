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
    
    
   
parser = argparse.ArgumentParser(description='Get small introns from UCSC mySQL database.')

parser.add_argument('-db', dest='database', default="hg19", type=str)
parser.add_argument('-table', dest='table', default="refGene", type=str)

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
        
        query = "SELECT * FROM {} WHERE (name2 LIKE 'MIR%' OR name2 LIKE 'LET%') AND exonCount = 1".format (table)
               
        cur.execute(query)
        
                
        isql = sql(cur)
                
        row = cur.fetchall()
                
        for r in row:
            
            chr, start, end, strand, name = isql.get (r,"chrom"), isql.geti (r,"txStart"), isql.geti (r,"txEnd"), isql.get (r, "strand"), isql.get (r, "name2")             
                
            if not (chr, start, end, strand) in dDone:
                
                print "\t".join ([chr, str(start), str(end), name, "-1", strand])
                        
                dDone[(chr, start, end, strand)] = 1
        
            
        print >>sys.stderr, "Found", len(dDone), "unique coordinates"
    
except mdb.Error, e:
  
    print "Error %d: %s" % (e.args[0], e.args[1])
    sys.exit(1)

finally:
    
    if con:
        con.close()

        
        
