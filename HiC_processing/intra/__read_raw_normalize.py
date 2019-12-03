#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 16:10:27 2018
Read chromosomal raw contact matrices from (Rao SSP 2014 )
(A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping)
normalize and correct for expected (O/E matrices)
for both intra and inter contact maps

@author: Naihui Zhou (nzhou@iastate.edu)
@updated 20190805
"""

import os
import numpy as np

chrs = [str(i) for i in range(1,23)]
chrs.extend('X')

class bin():
    def __init__(self, bin_left, chrm, step):
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
    
