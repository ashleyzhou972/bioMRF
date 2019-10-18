#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 16:10:27 2018
Read chromosomal raw contact matrices from (Rao SSP 2014 )
(A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping)
normalize and correct for expected (O/E matrices)
for both intra and inter contact maps

key info:
    cell type: IMR90
    resolution: 10kb
    MAPQ:>30
    normalization: KR
next step:
    map to ensembl genes
@author: Naihui Zhou (nzhou@iastate.edu)
@updated 20190805
"""

import os
import sys
import numpy as np

chrs = [str(i) for i in range(1,23)]
chrs.extend('X')
step = 10000

class bin():
    def __init__(self, bin_left, chrm):
        self.left = bin_left
        self.index = (int(bin_left)//step) 
        self.chrm = chrm


class interaction():
    def __init__(self, bin1, bin2):
        self.bin1 = bin1
        self.bin2 = bin2
        self.obs = 0
        self.norm = 0
    
    def read(self, interaction):
        self.obs = float(interaction)
    
    def normalize(self, norm_vec1, norm_vec2):
        """
        if normalization vector has NaN
        use 1 instead
        if normlaization vector is empty (chr9 for IMR90_10kb)
        use 1 instead
        """
        
        if len(norm_vec1)==0 or np.isnan(norm_vec1[self.bin1.index]):
            #print("norm1 empty\n")
            norm1 = 1
        else:
            norm1 = norm_vec1[self.bin1.index]
        if len(norm_vec2)==0 or np.isnan(norm_vec2[self.bin2.index]):
            #print("norm2 empty\n")
            norm2 = 1
        else:
            norm2 = norm_vec2[self.bin2.index]
        #updated 20190805    
        self.norm = float(self.obs)/(norm1*norm2)
        return(self.norm)
        
        
def bin_equal(bin1, bin2):
    if bin1.chrm==bin2.chrm:
        if bin1.left==bin2.left:
            return(True)
        else:
            return(False)
    else:
        return(False)



def read_raw_intra(rootfolder, resolution, chrm):
    """
    each bin is identified by chromosome and left side of bin
    """
    intra_list = list()
    raw_obs_file = os.path.join(rootfolder,'chr%s_%s.RAWobserved' % (chrm, resolution))
    
    with open(raw_obs_file, 'r') as f:
        for line in f:
            fields = line.strip().split()
            bin1 = bin(fields[0], chrm)
            bin2 = bin(fields[1], chrm)
            inter = interaction(bin1, bin2)
            inter.read(fields[2])
            intra_list.append(inter)                
    return(intra_list)
            
    
def read_raw_inter(rootfolder, resolution, chrm1, chrm2):
    """
    each bin is identified by chromosome and left side of bin
    """
    inter_list = list()
    if chrm1=='X' or (chrm2!='X' and int(chrm1)>=int(chrm2)):
        sys.stderr.write('Wrong chromosomal input for interchromosomal\n')
        sys.exit(1)
    inter_obs_file = os.path.join(rootfolder,'chr%s_%s_%s.RAWobserved' % (chrm1, chrm2, resolution))
    with open(inter_obs_file,'r') as f:
        for line in f:
            fields = line.strip().split()
            bin1 = bin(fields[0], chrm1)
            bin2 = bin(fields[1], chrm2)
            inter = interaction(bin1, bin2)
            inter.read(fields[2])
            inter_list.append(inter)
    return(inter_list)
   

    
def read_normalization_vector(rootfolder, norm_name, resolution, chrm):
    """
    norm_name is the name of the normalization method
    choice of KR, VC, SQRTVC
    """
    
    file = os.path.join(rootfolder,'chr%s_%s.%snorm' % (chrm, resolution, norm_name))
    norm_vec = np.fromfile(file,sep='\n')
    return(norm_vec)



def normalize(int_list, rootfolder, norm_name, resolution):
    chrm1 = int_list[0].bin1.chrm
    chrm2 = int_list[0].bin2.chrm
    norm_vec1 = read_normalization_vector(rootfolder, norm_name, resolution, chrm1)
    norm_vec2 = read_normalization_vector(rootfolder, norm_name, resolution, chrm2)
    for interaction in int_list:
        interaction.normalize(norm_vec1, norm_vec2)
    return(int_list)
    

"""
Below are legacy code, copied here, for importing purposes
original code can be found at /home/nzhou/hic/IMR90/work/read.py
"""

class coordinates():
    #build a class of genome coordinates
    def __init__(self, coords):
        #coords should be a tuple
        self.chrms = chrs 
        self.chrm = coords[0]
        self.start = coords[1]
        self.end = coords[2]
    def check(self):
        if type(self.start)==int and type(self.end)==int:
            if self.start<=self.end:
                if self.chrm in self.chrms:
                    return(True)
                else:
                    return(False)
            else:
                return(False)
        else:
            return(False)
            
            
def check_overlap(coords1, coords2):
    #coords1 should be locus coord (same size everywhere)
    #coords2 should be gene coord
    #Three states: disjoint, partial_overlap and complete overlap
    #If disjoint, return tuple(0)
    #If partial overlap, return tuple(1,percent_overlap), percent of coords1 that's overlapping with coords2
    #If complete coverage, return tuple(2,1): coords1 inside coords2 OR tuple(2,2) coords2 inside coords1
    if coords1.chrm==coords2.chrm:
        if coords1.end<coords2.start or coords1.start>coords2.end:#disjoint
            returnstate = 0
            returnnumber = None
        elif coords1.start<coords2.start and coords1.end>=coords2.start and coords1.end<=coords2.end :#partial
            returnstate = 1
            returnnumber = float(coords1.end-coords2.start)/float(coords1.end-coords1.start)
        elif coords1.start>coords2.start and coords1.start<=coords2.end and coords1.end>=coords2.end: #partial
            returnstate = 1
            returnnumber = float(coords2.end-coords1.start)/float(coords1.end-coords1.start) #percentage is with regard to coords1!
        elif coords1.start<=coords2.start and coords1.end>=coords2.end: #complete
            #coords1 covers coords2
            returnstate = 2
            returnnumber = float(coords2.end-coords2.start)/float(coords1.end-coords1.start)
        elif coords1.start>=coords2.start and coords1.end<=coords2.end: #complete
            #coords2 covers coords1
            returnstate = 2
            returnnumber = 1
        else:
            print('%s\t%s\n' % (coords1.start, coords1.end))
            print('%s\t%s\n' % (coords2.start, coords2.end))
    else:
        returnstate = 0
        returnnumber = None
    return(returnstate,returnnumber)




def convert_bin_to_coordinates(bin_object, step):
    coord = rrn.coordinates((bin_object.chrm, int(bin_object.left), int(bin_object.left)+step))
    return(coord)


def convert_coordinates_to_bin(coord, step):
    """
    coord must have
    coord.start % 500000==0
    """
    if coord.start % step==0:
        rbin = rrn.bin(str(coord.start), coord.chrm)
    else:
        rbin = None
        sys.stderr.write('Cannot convert coord to bin\n')
    return(rbin)
