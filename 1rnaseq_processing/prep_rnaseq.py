# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:44:30 2017
read ENCODE RNAseq quantification file
@author: nzhou
"""

import csv
import sys
import argparse
import glob
import os

def read_file(filename):
    countdict = dict()
    fpkmdict = dict()
    handle = open(filename,'r')
    handle.readline()
    for line in handle:
        fields = line.strip().split('\t')
        gene = fields[0]
        count = fields[7]
        fpkm = fields[6]
        countdict[gene]=float(count)
        fpkmdict[gene]=fpkm
    handle.close()
    return(countdict, fpkmdict)
    
def read_ensembl_mapping(mappingfile, coords):
    #coords should be a binary value
    #if True, return dictionary with ensembl gene id and genomic coordinates: (chrom, start, end)
    #if False, return dictionary with entrez gene id and ensembl gene id
    dictionary = dict()
    with open(mappingfile) as f:
        csvreader = csv.reader(f)
        next(csvreader)
        for line in csvreader:
            ensembl = line[0]
            chrom = line[1]
            start = line[2]
            end = line[3]
            entrez = line[4]
            if coords:
                dictionary[ensembl] = (chrom, start, end)
            else:
                dictionary[entrez] = ensembl
    return(dictionary)
                
        
        
def gene_id_cleaning(countdict, entrez_mapping, coords_mapping, filepath):
    #count UNIQUE gene IDs
    newdict = dict()
    deletedgene = set()
    for geneid in countdict.keys():
        #print(geneid)
        if geneid.startswith('EN'):
            newid = geneid.split('.')[0]
        else:
            try:
                #convert some entrez IDs to ensembl IDs
                newid = entrez_mapping[geneid]
            except KeyError:
                sys.stderr.write('NCBI gene ID %s not found for human.\n' % geneid)
                newid = None
        if newid!=None:
            try:
                coords_tuple = coords_mapping[newid]
            except KeyError:
                #sys.stderr.write('%s is not a protein-coding gene\n' % newid)
                deletedgene.add(newid)
                coords_tuple = None
            if coords_tuple!=None:
                if newid not in newdict.keys():
                    #value should be a tuple: (rnaseq_count, chrm, start, end)
                    newdict[newid] = (countdict[geneid],)+coords_tuple
                else:
                    #If there already exists a record for this id, happens from converting some entrez ids to ensembl id
                    #take the one with count that's larger
                    if float(newdict[newid][0])<float(countdict[geneid]):
                        newdict[newid] = (countdict[geneid],)+coords_tuple
    with open(filepath,'w') as f:
        f.write('%s\t%s\t%s\t%s\t%s\n' % ('Ensembl_ID','count','chromosome', 'start', 'end'))
        for gene in newdict:
            f.write('%s\t%s\t%s\t%s\t%s\n' % (gene, newdict[gene][0], newdict[gene][1], newdict[gene][2], newdict[gene][3])) 
    return(deletedgene)


    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description="Process ENCODE RNA-seq quantification")
    parser.add_argument("ensembl", type=os.path.abspath, help="File path to the Ensembl mapping file.")
    parser.add_argument("encode", type = os.path.abspath, help = "File path to the ENCODE gene quantification. Wildcard accepted")
    parser.add_argument("out", type=os.path.abspath, help = "Output folder")
    chrs = [str(i) for i in range(1,23)]
    chrs.extend('X')
    
    args = parser.parse_args()
    
    #1. Read ENSEMBL mapping file
    mappingfile = args.ensembl
    mapping_entrez = read_ensembl_mapping(mappingfile,False)
    mapping_coords = read_ensembl_mapping(mappingfile,True)
    #2.read raw gene quantification file
    for filename in glob.glob(args.encode):   
        encode_id = os.path.splitext(os.path.basename(filename))[0]
        print(encode_id)
        count, fpkm = read_file(filename)    
        print(args.out)
        outname = "rnaseq_%s.txt" % encode_id
        deleted = gene_id_cleaning(count, mapping_entrez, mapping_coords, os.path.join(args.out, outname))
