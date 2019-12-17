#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:41:39 2019

Map genes to TADs called by Arrowhead in Rao 2014.
TADs are intra-chromosomal structures!

@author: nzhou
@updated: 20190916
"""

import os
import sys
import argparse



class domain():
    def __init__(self, chrm1, start1, stop1, chrm2, start2, stop2, index):
        if chrm1==chrm2:
            self.chrm = chrm1
        else:
            sys.exit('Domain list contains inter-chromosomal domains\n')
        if start1==start2 and stop1==stop2:
            self.start = int(start1)
            self.stop = int(stop1)
        else:   
            sys.exit('Domain list contains x/y mismatch\n')
        self.index = index


    
def read_domainlist(domainlist_txt, chrm):
    domain_map = list()
    with open(domainlist_txt,'r') as dl:
        dl.readline()
        index = 0
        for line in dl:
            fields = line.strip().split('\t')
            chrm1 = fields[0]
            start1 = fields[1]
            stop1 = fields[2]
            chrm2 = fields[3]
            start2 = fields[4]
            stop2 = fields[5]
            # score = fields[7]
            if (chrm1==chrm and chrm2==chrm):
                index += 1
                dom = domain(chrm1, start1, stop1, chrm2, start2, stop2, index)
                domain_map.append(dom)
    return(domain_map)
            
    
    

def gene_mapping(rnaseq_file, domain_list, chrm):
    #chrm must agree with domain_list chrm!!!
    count = 0
    tad_genes = set()
    nontad_genes = set()
    with open(rnaseq_file,'r') as rna:
        rna.readline()
        for line in rna:
            fields = line.strip().split('\t')
            gene = fields[0]
            #print('gene: %s\n' % gene)
            gchrm = fields[2]
            if gchrm=='Y':
                #outhandle.write('%s\n' % gene) 
                continue
            if gchrm==chrm:
                count += 1
                gstart = int(fields[3])
                gend = int(fields[4])
                tad = False
                for dom in domain_list:
                    if gstart>=dom.start and gend<dom.stop:
                        tad_genes.add((gene, dom.index))
                        tad = True
                        break
                if not tad:
                    nontad_genes.add(gene)
        print("Total number of genes in chromosome %s is %d\n" % (chrm, count))
    return(tad_genes, nontad_genes)
                    
                   
            
            
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description="Read arrowhead TAD list and output whether genes are in TADs or not")
    
    parser.add_argument("domain_list_file", type=str, help = "File path to domain list")
    parser.add_argument("rnaseq_file", type = str, help = "File path to the RNAseq file with gene coordinates")
    parser.add_argument("outfolder", type=str, help="Path to output folder")
    
    chrs = [str(i) for i in range(1,23)]
    chrs.extend('X')
    
    parser.add_argument("--chr", dest = 'chrs', type = str, help = "Enter the chromosome, default is all.", default = chrs, choices=chrs)
    args = parser.parse_args()
    chrs = list(args.chrs)

    
    rnaseq_file = args.rnaseq_file
    outfolder = args.outfolder
    for chrm in chrs:
        domain_list=read_domainlist(args.domain_list_file,chrm)
        print("Number of TAD in chrm %s: %s" % (chrm, len(domain_list)))
        tad_genes, nontad_genes = gene_mapping(rnaseq_file, domain_list, chrm)
        print("Number of genes per TAD in chrm %s: %s\n" % (chrm, len(tad_genes)/len(domain_list)))
        
        with open(os.path.join(outfolder, 'chrm%s_TAD.txt' % chrm),'w') as tadout:
            for gene in tad_genes:
                tadout.write('%s\t%s\n' % (gene[0], gene[1]))
        with open(os.path.join(outfolder, 'chrm%s_nonTAD.txt' % chrm),'w') as nontadout:
            for gene in nontad_genes:
                nontadout.write('%s\n' % gene)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
