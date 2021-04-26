'''
    This script is for scoring the sequences based on
    four indicators.
    
    Written by Xiujia Yang, on July 2020
'''

import csv, os, glob, gzip, sys
import Levenshtein, argparse, subprocess, time
import pandas as pd
import seaborn as sns
import numpy as np
from Bio import SeqIO
from itertools import combinations
import matplotlib.pyplot as plt
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sklearn import linear_model
from math import sqrt
from numpy import mean
from Bio import pairwise2

def dist(seq1, seq2):
    '''
        This function is for comparing the distance
        of two zipped sequences
    '''
    dist = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != c2:
            dist += 1
    return dist

def cal_identity(allele , seq, seq_dict):
    '''
        This function accept the allele and discovered sequence
        and output the similarity to known sequences.
    '''
    if seq_dict.has_key(allele):
        if seq == seq_dict[allele]:
            iden_lst = []
            # Replace with the most prevalent length based on all allele corresponding to the same gene
            known_allele_lst = [ alt for alt in seq_dict.keys() if alt.split("*")[0] == allele.split("*")[0] and alt != allele ]
            # Replace with the most prevalent length based on all allele corresponding to the genes of the same family
            if known_allele_lst == []:
                known_allele_lst = [ alt for alt in seq_dict.keys() if alt.split("-")[0] == allele.split("-")[0] ]
            # Replace with the most prevalent length based on all allele
            if known_allele_lst == []:
                known_allele_lst = seq_dict.keys()
            for alt in known_allele_lst:
                seq1 = seq
                seq2 = seq_dict[alt]
                #seq1, seq2 = pairwise2.align.globalxx(seq1, seq2)[0][:2]
                seq1, seq2 = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -1)[0][:2]
                iden = 1-float(dist(seq1, seq2))/len(seq1)
                iden_lst.append(iden)
            return mean(iden_lst)     
        else:
            seq1 = seq
            seq2 = seq_dict[allele]
            seq1, seq2 = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -1)[0][:2]
            iden = 1-float(dist(seq1, seq2))/len(seq1)
            return iden
    else:
        iden_lst = []
        # Replace with the most prevalent length based on all allele corresponding to the same gene
        known_allele_lst = [ alt for alt in seq_dict.keys() if alt.split("*")[0] == allele.split("*")[0] ]
        # Replace with the most prevalent length based on all allele corresponding to the genes of the same family
        if known_allele_lst == []:
            known_allele_lst = [ alt for alt in seq_dict.keys() if alt.split("-")[0] == allele.split("-")[0] ]
        # Replace with the most prevalent length based on all allele
        if known_allele_lst == []:
            known_allele_lst = seq_dict.keys()
    
        for alt in known_allele_lst:
            seq1 = seq
            seq2 = seq_dict[alt]
            #seq1, seq2 = pairwise2.align.globalxx(seq1, seq2)[0][:2]
            seq1, seq2 = pairwise2.align.globalms(seq1, seq2, 2, -1, -5, -1)[0][:2]
            iden = 1-float(dist(seq1, seq2))/len(seq1)
            iden_lst.append(iden)
        return mean(iden_lst)

def main():
    ## obtain the reference sequence and store it in dictionary
    seq_dict = {}
    for rec in SeqIO.parse(refFl, "fasta"):
        if "partial" not in rec.description:
            allele = rec.id.split("|")[1]
            seq_dict.update({allele: str(rec.seq).upper()})
    
    ## import the input material files
    df = pd.read_csv(seqFl, sep='\t', index_col=0)
    
    ## Calculate the constant values to serve as denominators
    species = seqFl[2:].split("_")[0]
    chain = seqFl.split("_")[1]
    if species == "rhesus" and chain == "igk":
        RT = df[df.IMGT_leader==4].groupby("seq_ref")["readFreq"].sum().quantile(q=0.50)
        CT = df[df.IMGT_leader==4].groupby("seq_ref")["cloneFreq"].sum().quantile(q=0.50)
    else:
        RT = df[df.IMGT_leader==0].groupby("seq_ref")["readFreq"].sum().quantile(q=0.50)
        CT = df[df.IMGT_leader==0].groupby("seq_ref")["cloneFreq"].sum().quantile(q=0.50)
    
    
    ## Iterate over each sequence to calculate its score
    out = csv.writer(open(outFl, "wb"), delimiter="\t")
    title = ["Allele", "Sequence", "nDonors", "nClones", "nReads", "Identity", "Sd", "Sc", "Sr", "Ss", "Stotal", "Flag", "leaderStart"]
    out.writerow(title)
    for seq_allele, group in df.groupby("seq_ref"):
        seq, allele = seq_allele.split("_")
        leader_start = group.leaderStart.tolist()[0]
        nReads = group["readFreq"].sum()
        nClones = group["cloneFreq"].sum()
        nDonors = group.shape[0]
        Sr = 20 * sqrt(nReads/float(RT)); Sr = 20 if Sr > 20 else Sr
        Sc = 20 * sqrt(nClones/float(CT)); Sc = 20 if Sc > 20 else Sc
        iden = cal_identity(allele, seq[leader_start:], seq_dict)
        Ss = 30 * iden; Ss = 30 if Ss > 30 else Ss
        Sd = 30 * (1-1/(float(nDonors+1))); Sd = 30 if Sd > 30 else Sd
        Stotal = Sr + Sc + Sd + Ss
        Flag = group.IMGT_leader.tolist()[0]
        outLine = [allele, seq, nDonors, nClones, nReads, iden, Sd, Sc, Sr, Ss, Stotal, Flag, leader_start]
        out.writerow(outLine)
    
if __name__ == "__main__":
    seqFl = sys.argv[1]  # discovered sequence list
    refFl = sys.argv[2]  # reference leader or utr sequences
    outFl = sys.argv[3]  # the name of output score file
    ## Run main function
    main()
