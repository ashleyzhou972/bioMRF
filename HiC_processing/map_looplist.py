#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 12:43:06 2019
Map loop list to genes
loop list downloaded from Rao et al. 2014
@author: nzhou
"""

"""
Downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
filename: GSE63525_IMR90_HiCCUPS_looplist.txt.gz
key info:
    cell type: IMR90
    resolution: 10kb
"""
import os
os.chdir('/home/nzhou/hic/rao2014')
import read_raw_normalize as rrn
chrms = [str(i) for i in range(1,24)]
chrms.append('X')

def read_and_write_looplist(looplist_file, folder):
    """
    See https://github.com/aidenlab/juicer/wiki/HiCCUPS#loop-list-content
    for file format
    Write looplist file for each chromosome
    Include coordinates and centroid coordinates
    """
    chrms = [str(i) for i in range(1,24)]
    chrms.append('X')
    handles = dict()
    with open(looplist_file,'r') as ll:
        header = ll.readline()
        header_fields = header.strip().split('\t')
        sub_header = [header_fields[i] for i in [0,1,2,3,4,5,7,17,18,19]]
    for chrm in chrms:
        handle =open(os.path.join(folder, 'looplist_chrm%s.txt' % chrm), 'w')
        handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % tuple(sub_header))
        handles[chrm] = handle
    with open(looplist_file,'r') as ll:
        header = ll.readline()
        for line in ll:
            fields = line.strip().split('\t')
            sub_fields = [fields[i] for i in [0,1,2,3,4,5,7,17,18,19]]
            chrm = fields[0]
            handles[chrm].write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % tuple(sub_fields))
    for chrm in chrms:
        handles[chrm].close()
  
      
def read_looplist(folder, chrm_in):
    """
    Read chromosome-specific loop list with 10 fields instead of the original 20
    file name is 'looplist_chrm%s.txt
    @param chrms a list of chromosomes of interest. chromosome number need to be char'
    """
    ret = list()
    inhandle = open(os.path.join(folder, 'looplist_chrm%s.txt' % chrm_in), 'r')
    next(inhandle)
    for line in inhandle:
        fields = line.strip().split('\t')
        xcoords = rrn.coordinates((fields[0], int(fields[1]), int(fields[2])))
        ycoords = rrn.coordinates((fields[3], int(fields[4]), int(fields[5])))
        paircoords = (xcoords, ycoords)
        ret.append(paircoords)
    inhandle.close()
    return(ret)

    
def gene_mapping_looplist(rnaseq_file, looplist_list, chrm_in):       
    """
    for each gene, map it to a bin, save in a dictionary
    copied from ./map_ensembl_genes.py
    modified on 20190409
    @param rnaseq_file file that contains rnaseq data with gene name, and location information
    @param looplist_list a list of loops represented by tuple of two loci, each represented with rrn.coordinates class
    @param chrms a list of chromosomes of interest. chromosome number need to be char
    chrms need to be identical with the chrms that generated looplist_list
    """
    gene_dict = dict() #dict containing the most
    with open(rnaseq_file,'r') as rna:
        rna.readline()
        for line in rna:
            fields = line.strip().split('\t')
            gene = fields[0]
            #print('gene: %s\n' % gene)
            chrm = fields[2]
            if chrm=='Y':
                continue
            if chrm!=chrm_in:
                continue
            start = int(fields[3])
            end = int(fields[4])
            genecoords = rrn.coordinates((chrm,start,end))
            if genecoords in gene_dict.keys():
                print("duplicate gene record%s\t" % gene)
            gene_dict[genecoords] = gene
    results = list()
    index = 0
    for loop in looplist_list: 
        #each loop has two loci
        #i = 0 for upstream
        #i = 1 for downstream
        locus0 = loop[0]
        locus0_result =list() 
        for genecoords in gene_dict:
            overlap_state, overlap_number = rrn.check_overlap(locus0, genecoords)
            if overlap_state!=0:
                locus0_result.append((gene_dict[genecoords], overlap_state, overlap_number))
        locus1 = loop[1]
        locus1_result =list() 
        for genecoords in gene_dict:
            overlap_state, overlap_number = rrn.check_overlap(locus1,genecoords)
            if overlap_state!=0:
                locus1_result.append((gene_dict[genecoords], overlap_state, overlap_number))
        loop_result = (index, locus0_result, locus1_result)
        results.append(loop_result)
        index+=1
    return(results)
    
def get_higher_overlap(locus_gene_result, threshold):
    #input should be list of tuples
    #Each tuple should be of three elements: genename, overlap_state(0,1,2) and overlap_number
    #return one out of all the tuples
    #if there is no tuple, then return an empty list
    #if there is one tuple then return that one
    
    #overlap state 1: partial 2:complete
    #threshold for overlap: if larger than threshold count as overlap, otherwise not
    if len(locus_gene_result)==0:
        return(())
    elif len(locus_gene_result)==1:
        gene_state = locus_gene_result[0][1]
        gene_overlap = locus_gene_result[0][2]
        if gene_state==2: #complete  overlap
            return(locus_gene_result[0])
        else:
            if gene_overlap<=threshold: #overlap less than threshold
                return(())
            else:
                return(locus_gene_result[0])
    else:
        max_state = 1
        max_overlap = threshold
        max_gene = ()
        for gene in locus_gene_result:
            gene_state = gene[1]
            gene_overlap = gene[2]
            if gene_state>max_state:
                max_gene = gene
                max_state = gene_state
            elif gene_state==max_state:
                if gene_overlap > max_overlap:
                    max_gene = gene
                    max_overlap = gene_overlap
        return(max_gene)
                
    
    
def write_edge_list(loopmap, outfolder, chrm, threshold):
    nomap_count = 0
    indices = set()
    note = open(os.path.join(outfolder, 'edge_list_with_overlap_chrm%s.txt' % chrm), 'w')
    with open(os.path.join(outfolder, 'edge_list_chrm%s.txt' % chrm),'w') as of:
        for index in range(0,len(loopmap)):
            locus0 = loopmap[index][1]
            locus1 = loopmap[index][2]
            final0 = get_higher_overlap(locus0, threshold)
            final1 = get_higher_overlap(locus1, threshold)
            if final0 and final1:
                #if they are both not empty
                of.write('%s\t%s\n' % (final0[0], final1[0])) #only write the gene names
                indices.add(index)
            else:
                nomap_count += 1
                if not final0: 
                    final0 = ('','','')
                if not final1:
                    final1 = ('','','')
            note.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (index, final0[0], final0[1], final0[2],
                                         final1[0], final1[1], final1[2]))
    note.close()
    return(indices, nomap_count)
            
            
    
    
if __name__=='__main__':
    #looplist_file ='/home/nzhou/hic/rao2014/GSE63525_IMR90_HiCCUPS_looplist.txt'
    #read_and_write_looplist(looplist_file, folder)
    folder='/home/nzhou/hic/rao2014/res10k'
    subfolder = os.path.join(folder, 'thres10')
    rnaseq_file ='./rnaseq_corrected.txt'
    for chrm in chrms:
        looplist = read_looplist(os.path.join(folder, 'looplists'), chrm)
        loop_maps= gene_mapping_looplist(rnaseq_file, looplist, chrm)
        indices, nomap_count = write_edge_list(loop_maps, subfolder, chrm, 0.1)
        with open(os.path.join(subfolder,'matched_loops_indices_chrm%s.txt' % chrm),'w') as ind_handle:
            for i in indices:
                ind_handle.write('%s\n' % i)
        
    
            