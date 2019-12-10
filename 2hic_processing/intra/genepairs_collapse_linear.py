#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:41:44 2019
get edge lists for linear neighborhood networks
#updated 201906
@author: nzhou
"""


import os
import numpy as np
import argparse

def read_and_sort_genes(rnaseq_file, chrs):
    master_dict1 = dict()
    master_dict2 = dict()
    for chrm in chrs:
        master_dict1[chrm] = dict()
        master_dict2[chrm] = dict()
    with open(rnaseq_file,'r') as rf:
        rf.readline()
        for line in rf:
            fields = line.strip().split('\t')
            gene = fields[0]
            #print('gene: %s\n' % gene)
            chrm = fields[2]
            if chrm=='Y':
                #outhandle.write('%s\n' % gene) 
                continue
            start = int(fields[3])
            end = int(fields[4])
            master_dict1[chrm][gene] = (start, end)
            master_dict2[chrm][gene] = start
    sorted_dict = dict()
    for chrm in chrs:
        sorted_dict[chrm] = [(k, master_dict1[chrm][k]) for k in sorted(master_dict2[chrm], key = master_dict2[chrm].get, reverse=False)]
    return(sorted_dict)
    
def get_gene_pairs(sorted_dict, outfolder, cutoff):
    for chrm in sorted_dict:
        outhandle = open(os.path.join(outfolder, "linear_chr%s.genepairs" % chrm),'w')
        sortedlist = sorted_dict[chrm]
        for i in range(0, len(sortedlist)-1):
            end = sortedlist[i][1][1]
            start = sortedlist[i+1][1][0]
            if start-end < cutoff:
                outhandle.write("%s\t%s\n" % (sortedlist[i][0], sortedlist[i+1][0]))
        outhandle.close()
            
def analyze_gene_pairs(sorted_dict, chrm):
    sortedlist = sorted_dict[chrm]
    dists = list()
    for i in range(0, len(sortedlist)-1):
        end = sortedlist[i][1][1]
        start = sortedlist[i+1][1][0]
        dist = start - end
        dists.append(dist)
    dists_arr = np.array(dists)
    print(np.median(dists_arr))
    return(dists_arr)
    
            
    
        
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description="Output gene pairs for linear gene network")
    parser.add_argument("root", type=os.path.abspath, help="Root folder")
    parser.add_argument("coords_file", type = os.path.abspath, help = "File path to the processed RNA-seq file with gene coordinates")
    parser.add_argument("cutoff", type=int, help="Genomic distance cutoff for neighborhood inference. e.g. A cutoff of 10k means genes more than 10k bp apart are not neighbors")
    
    chrs = [str(i) for i in range(1,23)]
    chrs.extend('X')
    
    args = parser.parse_args()
    
    rnaseq_file = args.coords_file

    sorted_dict = read_and_sort_genes(rnaseq_file, chrs)

    get_gene_pairs(sorted_dict, args.root, args.cutoff)
