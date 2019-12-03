#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:17:20 2019

Map gene names from a gaf file 

@author: nzhou
@updated 20190828
"""

import os
import Bio.UniProt.GOA as GOA


def readGaf(gaf):
    prot_set = set()
    ingafhandle = open(gaf,'r', encoding = 'utf-8-sig')
    for rec in GOA.gafiterator(ingafhandle):
        if not rec['Qualifier'][0].startswith('NOT'):
            if rec['DB']=='UniProtKB':
                ac = rec['DB_Object_ID']
                prot_set.add(ac)
        
    ingafhandle.close()
    return(prot_set)



if __name__=="__main__":
    folder = '/home/nzhou/hic/rao2014/IMR90_10kb/mitochondria/'
    mito = os.path.join(folder, 'gaf_0005739_mito.gaf')
    prots = readGaf(mito)
    with open(os.path.join(folder, 'mito_genes_uniprot.txt'),'w') as out:
        for prot in prots:
            out.write("%s\n" % prot)
    ##Use ID mapping on UniProt to map to Ensembl