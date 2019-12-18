#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 16:12:59 2019

@author: Naihui Zhou (nzhou@iastate.edu)
@last updated 20191217
"""
import heapq
import os
import gzip
from Ontology.IO import OboIO
import Bio.UniProt.GOA as GOA
import argparse

"""
 First write ancestor files using obo
 Note that even though we are concerned with BPO only, 
 a lot of MFO terms have BPO terms as ancestors 
 so we are including MFO ancestor files as well

 Edges include 'part_of', 'is_a' and 
 AND 'positively_regulates', 'negatively_regulates' and 'regulates'!!!
"""

def go_ontology_split(ontology):
    """
    Split a GO obo object into three ontologies
    """
    mfo_terms = set({})
    bpo_terms = set({})
    cco_terms = set({})
    for node in ontology.get_ids(): # loop over node IDs and alt_id's
        if ontology.namespace[node] == "molecular_function":
            mfo_terms.add(node)
        elif ontology.namespace[node] == "biological_process":
            bpo_terms.add(node)
        elif ontology.namespace[node] == "cellular_component":
            cco_terms.add(node)
        else:
            raise(ValueError,"%s has no namespace" % node)
    return (mfo_terms, bpo_terms, cco_terms)

def go_ontology_ancestors_split_write(obo_path):
    """
    Input: an OBO file
    Output: 3 files with ancestors
    """
    
    obo_bpo_out = open("%s_ancestors_bpo.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_cco_out = open("%s_ancestors_cco.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_mfo_out = open("%s_ancestors_mfo.txt" % (os.path.splitext(obo_path)[0]),"w")
    obo_parser = OboIO.OboReader(open(obo_path))
    go = obo_parser.read()
    mfo_terms, bpo_terms, cco_terms = go_ontology_split(go)
    for term in mfo_terms:
        ancestors = go.get_ancestors(term)
        obo_mfo_out.write("%s\t%s\n" % (term,",".join(ancestors)))
    for term in bpo_terms:
        ancestors = go.get_ancestors(term)
        obo_bpo_out.write("%s\t%s\n" % (term,",".join(ancestors)))
    for term in cco_terms:
        ancestors = go.get_ancestors(term)
        obo_cco_out.write("%s\t%s\n" % (term,",".join(ancestors)))

    obo_mfo_out.close()
    obo_bpo_out.close()
    obo_cco_out.close()
    return([len(bpo_terms),len(cco_terms),len(mfo_terms)])


"""
 Read GOA annotations, propogate through ancestors
 Two different strategies: propogate leaf MFO terms to ancestor BPO terms
  vs. not propogate leafy MFO terms to ancestor BPO terms
 Notice the depth of the BPO terms
 
 Accepted evidence codes include: 
     'EXP','IDA','IPI','IMP','IGI','IEP','TAS','IC'
     See http://geneontology.org/docs/guide-go-evidence-codes/
""" 


def readAncestor(ancestorFile):
    ancestorDict = dict()
    with open(ancestorFile,'r') as af:
        for line in af:
            goterm = line.strip().split('\t')[0]
            try:
                ancestors = set(line.strip().split('\t')[1].split(','))
            except IndexError:
                ancestors = set()
            ancestorDict[goterm] = ancestors
    return(ancestorDict)
 

## Below function write leaf annotations to a tab file
def read_gaf_write_tab(gaf_file, include_mfo, outfile):
    Evidence = {'Evidence': set(['EXP','IDA','IPI','IMP','IGI','IEP','TAS','IC'])}
    if include_mfo:
        Aspect = {'Aspect':set(['P','F'])}
    else:
        Aspect = {'Aspect':set(['P'])}
    Evidence = {'Evidence': set(['EXP','IDA','IPI','IMP','IGI','IEP','TAS','IC'])}
    outhandle = open(outfile, 'w')
    ingafhandle = open(gaf_file,'r')
    counter = 0
    for rec in GOA.gafiterator(ingafhandle):
        if GOA.record_has(rec, Aspect):
            if GOA.record_has(rec, Evidence):
                prot = rec['DB_Object_ID']
                go = rec['GO_ID']
                outhandle.write("%s\t%s\n" % (prot, go))
                counter +=1
    ingafhandle.close()
    outhandle.close()
    return(counter)

#three root terms

#returns set of direct children of GO:0008150 Biological Process
def read_children(children_tab):
    direct_children = set()
    with open(children_tab, 'r') as tab:
        for line in tab:
            go = line.split()[0]
            direct_children.add(go)
    return(direct_children)
    
    
    
## returns true if the given go-term is not root,
## and belongs to the given aspect 
    
def check_go(go_term, go_object, roots, aspect):
    aspectDict = {"C":"cellular_component", "F": "molecular_function","P":"biological_process"}
    if go_term in roots:
        return(False)
    else:
        try:
            namespace = go_object.namespace[go_term]
        except KeyError:
            print("GO term %s not found in obo\n" % go_term)
            return(False)
        if namespace==aspectDict[aspect]:
            return(True)
        else:
            return(False)

# Below function returns 
# a dictionary with GO-term as key
# and set of prots as value
def read_tab_propogate(tab_file, ancestor_files, go_object, aspect, root_terms):
    outDict = dict()
    anceDict = dict()
    for anc_file in ancestor_files:
        anceDict = {**readAncestor(anc_file), **anceDict}
    with open(tab_file, 'r') as tab:
        for line in tab:
            prot, go = line.strip().split('\t')
            try:
                ancestors = anceDict[go]
            except KeyError:
                print("GO term %s is not found in OBO\n" % go)
                continue
            ancestors.add(go)  #add itself
            for term in ancestors:
                if check_go(term, go_object, root_terms, aspect):
                    if term not in outDict.keys():
                        outDict[term] = set()
                    outDict[term].add(prot)
    return(outDict)
            
    

def read_tab_no_propogate(tab_file, ancestor_file, go_object, aspect, root_terms):
    outDict = dict()
    anceDict = readAncestor(ancestor_file)
    print(len(anceDict))
    with open(tab_file, 'r') as tab:
        for line in tab:
            prot, go = line.strip().split('\t')
            term = go
            if check_go(term, go_object, root_terms, aspect):
                if term not in outDict.keys():
                    outDict[term] = set()
                outDict[term].add(prot)
    return(outDict)
## Below filters the Uniprot proteins (gene products) according to our
## gene list in ensembl
## returns a set of not mapped ensembl genes
## and a dictionary with keys being ensembl genes, 
## and values being sets of prots in uniprot accession
def filter_genes(mapping_file, ensembl_genes):
    genes = set()
    with open(ensembl_genes,'r') as eg:
        for line in eg:
            gene = line.strip()
            genes.add(gene)
    print("Total number of genes queried: %d\n" % len(genes))
    mapping_dict = dict()
    with gzip.open(mapping_file, 'rt') as mf:
        mf.readline()
        for line in mf:
            ensembl_id = line.strip().split()[0]
            uniprot_id = line.strip().split()[3]
            if ensembl_id in genes:
                #remove from set of genes to see unmapped genes
                genes.remove(ensembl_id)
                if uniprot_id not in mapping_dict.keys():
                    mapping_dict[uniprot_id] = set()
                mapping_dict[uniprot_id].add(ensembl_id)
    print("Number of genes not mapped: %d\n" % len(genes))
    return(genes, mapping_dict)
      

##Get the top BPO terms that annotated the most proteins
def get_top_terms(goa_dict, N):
    mylist = list()
    for key in goa_dict:
        heapq.heappush(mylist, (len(goa_dict[key]), key))
    return(heapq.nlargest(N, mylist))
    
    
#output the list of ensembl genes associated with a GO term
def get_ensembl_genes(go_term_list, goa_dict, ensembl_mapped_dict):
    ret_dict = dict()
    for go_term in go_term_list:
        ensembl = set()
        for prot in goa_dict[go_term]:
            try:
                genes = ensembl_mapped_dict[prot]
                ensembl = ensembl.union(genes)
            except KeyError:
                #print("Protein %s not mapped to genes\n" % prot)
                #This should be a common occurence so can suppress warning print
                pass
        ret_dict[go_term] = ensembl
    return(ret_dict)
    
## translate GO terms to its names
def go_term_translator(go_terms, go_object):
    names = dict()
    for term in go_terms:
        name = go_object.get_term(term).name
        names[term] = name
    return(names)


if __name__=="__main__":
    
    ##argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--propagate', dest='propogate', action='store_true')
    parser.add_argument('--no-propagate', dest='propogate', action='store_false')
    parser.add_argument("N", help="choose the top N most used GO terms", type=int)
    parser.add_argument("datafolder",type = str, help="Path to folder where input and output data are saved\n")
    parser.add_argument("genename_folder",type = str, help="Path to gene name folder\n")
    args = parser.parse_args()
    
    
    exclude_direct_children_of_root = False

    datafolder = args.datafolder
    obo_file = os.path.join(datafolder, "go.obo")
    obo_parser = OboIO.OboReader(open(obo_file,'r'))
    go_object = obo_parser.read()
    
    #c_bpo, c_cco, c_mfo = go_ontology_ancestors_split_write(obo_file)
    
    #gaf_file = os.path.join(datafolder, "goa_human.gaf")
    #c_nomfo = read_gaf_write_tab(gaf_file, False, os.path.join(datafolder, "goa_human_bpo.tab"))
    #c_mfo = read_gaf_write_tab(gaf_file, True, os.path.join(datafolder, "goa_human_bpo_mfo.tab"))


    ##map Uniprot proteins to ensembl genes
    ensembl_file = os.path.join(args.genename_folder, 'all_gene_names.txt')
    ## We use the mapping file from the Ensembl 90 release 
    ## Because some of our genes are obsolete in the newer version!!
    mapping_file = os.path.join(datafolder, 'Ensembl_uniprot.tsv.gz')
    notmapped_set, mapped_dict = filter_genes(mapping_file, ensembl_file)
    
    
    ## List of terms excluded (due to being too general)    
    root_terms = set(['GO:0008150', 'GO:0005575','GO:0003674'])
    if exclude_direct_children_of_root:
        direct_children = read_children(os.path.join(datafolder, 'direct_children_of_bpo.tab'))
        to_exclude = root_terms.union(direct_children)
    else:
        to_exclude = root_terms
    

    ## bpo no prop
    if args.propogate:
        tab_file_mfo_bpo = os.path.join(datafolder, "goa_human_bpo_mfo.tab")
        ancestor_file_bpo = os.path.join(datafolder, "go_ancestors_bpo.txt")
        ancestor_file_mfo = os.path.join(datafolder, "go_ancestors_mfo.txt")
        anno = read_tab_propogate(tab_file_mfo_bpo, [ancestor_file_bpo, ancestor_file_mfo], go_object, 'P', to_exclude)
    else:
        tab_file_bpo = os.path.join(datafolder, "goa_human_bpo.tab")
        ancestor_file_bpo = os.path.join(datafolder, "go_ancestors_bpo.txt")
        anno = read_tab_no_propogate(tab_file_bpo, ancestor_file_bpo, go_object, 'P', to_exclude)
    ## output the n terms annotating the most proteins
    N = args.N
    query_terms = [x[1] for x in get_top_terms(anno, N)]
    query_names = go_term_translator(query_terms, go_object)
    
    
    ## Output the ensembl genes associated with each of the N GO terms
    ensembl_list_dir = os.path.join(datafolder, "ensembl_list")
    if not os.path.exists(ensembl_list_dir):
        os.makedirs(ensembl_list_dir)
    output = get_ensembl_genes(query_terms, anno, mapped_dict)
    for goterm in output:
        with open(os.path.join(datafolder, "ensembl_list", "%s_%s_genes_ensembl.txt" % (goterm, args.propogate)), 'w') as out:
            for ensembl in output[goterm]:
                out.write("%s\n" % ensembl)
                
                
    ## write an outfile of all the top N GO terms and their prot count and gene

    with open(os.path.join(datafolder, "Top%d_GO_terms_counts_%s" % (N, args.propogate)),"w") as out_count:
        for term in get_top_terms(anno, N):
            out_count.write("%s\t%s\t%s\t%s\n" % (term[1], query_names[term[1]], term[0], len(output[term[1]])))
    
