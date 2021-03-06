#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 14:21:47 2018

read ensembl gene id coordinates and map them to contact matrix bins
This is regardless of HiC data
Only hyperparameter is the size of the bin

@author: Naihui Zhou nzhou@iastate.edu
@update: 20191122

"""

import os
import __read_raw_normalize as rrn
import argparse


def compare_chrms(chrm1, chrm2):
    """
    return true if chrm1 < chrm2
    """
    if chrm1=='X':
        return(False)
    else:
        if chrm2=='X':
            return(True)
        else:
            if int(chrm1)<int(chrm2):
                return(True)
            else:
                return(False)

def mathematical_bin_overlap(left, right, step, threshold):
    """
    instead of going through all possible bins, 
        calculate the quotient and remainder of the left and right 
        coordinates of genes
    """
    div_left = divmod(left, step)
    div_right = divmod(right, step)
    if div_left[1]<=(1-threshold)*step:
        index_left = div_left[0]
    else:
        index_left = div_left[0]+1
    if div_right[1]>threshold*step:
        index_right = div_right[0]
    else:
        index_right = div_right[0]-1
    return([x for x in range(index_left*step, (index_right+1)*step,step)])
        
    
        
    
    
def gene_mapping_nodup(rnaseq_mapping_file,  threshold):
    """
    updated 20190421
    for each gene, map it to all the bins that have at least 
        "threshold"% overlap. save in dictionary
    updated 20190424
    Do not go over all bin lefts
    Simply do divison and remainders
    """
    gene_dict = dict()
    #outhandle = open(outfile_Ychr,'w')
    with open(rnaseq_mapping_file,'r') as rna:
        rna.readline()
        for line in rna:
            fields = line.strip().split('\t')
            gene = fields[0]
            #print('gene: %s\n' % gene)
            chrm = fields[2]
            if chrm=='Y':
                #outhandle.write('%s\n' % gene) 
                continue
            start = int(fields[3])
            end = int(fields[4])
            bins = mathematical_bin_overlap(start, end, step, threshold)
            result_bins = set(bins)
            gene_dict[gene] = result_bins
    #outhandle.close()
    return(gene_dict)

def read_gene_chromosome(rnaseq_mapping):
    chrm_dict = dict()
    with open(rnaseq_mapping,'r') as r:
        for line in r:
            fields = line.strip().split('\t')
            gene = fields[0]
            chrm = fields[2]
            chrm_dict[gene]=chrm
    return(chrm_dict)
    
    
def reverse_dict(gene_dict, chrm_dict, chrm):
    #updated 20190424: each gene maps to multiple bins
    bin_dict = dict()
    for gene in gene_dict:
        chrm_gene = chrm_dict[gene]
        if chrm_gene==chrm:
            for bins in gene_dict[gene]:
                if str(bins) not in bin_dict.keys():
                    bin_dict[str(bins)] = set()
                bin_dict[str(bins)].add(gene)
        else:
            continue
    return(bin_dict)

def convert_resolution_to_step(reso):
    step = int(reso[:(reso.find('kb'))])*1000
    return(step)

def inter_filter_normalize_and_write_raw(rootfolder, norm_name, resolution, chrm1, chrm2, bin_dict1, bin_dict2):
    """
    bin_dict must match chrm!!!
    each bin is identified by chromosome and left side of bin
    """
    outfile = open(os.path.join(rootfolder, 'norm','inter_chr%s_chr%s_%s.norm' % (chrm1, chrm2, resolution)),'w')
    outfile.write('%s\t%s\t%s\n' % ('gene1', 'gene2', 'inter_norm'))
    norm_vec1 = rrn.read_normalization_vector(os.path.join(rootfolder, 'raw'), norm_name, resolution, chrm1)
    norm_vec2 = rrn.read_normalization_vector(os.path.join(rootfolder, 'raw'), norm_name, resolution, chrm2)
    raw_obs_file = os.path.join(rootfolder,'raw', 'chr%s_%s_%s.RAWobserved' % (chrm1, chrm2, resolution))
    with open(raw_obs_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            left1 = fields[0]
            bin1 = rrn.bin(left1, chrm1, step)
            left2 = fields[1]
            bin2 = rrn.bin(left2, chrm2, step)
            if left1 in bin_dict1.keys() and left2 in bin_dict2.keys():
                inter = rrn.interaction(bin1, bin2)
                inter.read(fields[2])
                inter_norm = inter.normalize(norm_vec1, norm_vec2)
                gene1_set = bin_dict1[left1]
                gene2_set = bin_dict2[left2]
                for gene1 in gene1_set:
                    for gene2 in gene2_set:
                        outfile.write('%s\t%s\t%s\n' % (gene1, gene2, inter_norm))
    outfile.close()
    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Read raw HiC contact matrix and normalize")
    parser.add_argument("root", type=str, help="Root folder")
    parser.add_argument("coords_file", type = str, help = "File path to the Ensembl file with gene coordinates")
    parser.add_argument("--resolution", dest = 'reso', type=str, help = "Enter resolution, e.g. '10kb' or '250kb'.", default = "10kb")
    parser.add_argument("--overlap", dest = 'overlap', type=float, help="percentage of overlap between a gene and a locus", default = 0.1)
    parser.add_argument("--norm", dest = 'norm_name', help = "Normalization method", default = "KR", choices=["KR", "VC", "SQRTVC"]) 
    chrs = [str(i) for i in range(1,23)]
    chrs.extend('X')
    
    parser.add_argument("--chrA", dest = 'chrsA', type = str, help = "Enter the first chromosome in the pair, default is all.", default = chrs, choices=chrs)
    parser.add_argument("--chrB", dest = 'chrsB', type = str, help = "Enter the second chromosome in the pair, default is all.", default = chrs, choices=chrs)
    args = parser.parse_args()
    chrsA = list(args.chrsA)
    chrsB = list(args.chrsB)

    step = convert_resolution_to_step(args.reso)
    
    rnaseq_mapping_file = args.coords_file
    newtt = gene_mapping_nodup(rnaseq_mapping_file, args.overlap)
    chrm_dict = read_gene_chromosome(rnaseq_mapping_file)
    bin_dict_dict = dict()
    for chrm1 in chrsA:
        bin_dict1 = reverse_dict(newtt, chrm_dict, chrm1)
        bin_dict_dict[chrm1] = bin_dict1
        # Create dictionary for bin_dicts to avoid reading them repeatedly
        for chrm2 in chrsB:
            if compare_chrms(chrm1, chrm2):
                print("Processing Chromosome pair %s and %s\n" % (chrm1, chrm2))
                try:
                    bin_dict2 = bin_dict_dict[chrm2]
                except KeyError:
                    print("bin_dict for %s do not already exist\n" % chrm2)
                    bin_dict2 = reverse_dict(newtt, chrm_dict, chrm2)
                    bin_dict_dict[chrm2] = bin_dict2
                inter_filter_normalize_and_write_raw(args.root, args.norm_name, args.reso, chrm1, chrm2, bin_dict1, bin_dict2)
        
    

