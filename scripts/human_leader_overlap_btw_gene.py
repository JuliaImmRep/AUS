import csv, os, glob, gzip, sys
import pandas as pd
import seaborn as sns
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from sklearn import linear_model
import re
from math import sqrt, ceil
from collections import Counter
from scipy.stats import ttest_ind

def orphan_rename(name):
    '''
        rename the orphan gene to enable valid file output
    '''
    if "/OR" in name:
        return "-".join(name.split("/"))
    else:
        return name


def main():
    for chain in ["igh", "igk", "igl"]:
        ## Obtain the gene usage table with row being genes and columns being samples
        df_gu = pd.read_csv("data/human.%s.gene.expr.txt"%chain, sep="\t", index_col=0)
        df_gu = df_gu.stack().reset_index().rename(columns={'level_0':'gene','level_1':'sample', 0:'gene usage'})
        
        ## Obtain the filtered upstream sequences and extract from each sample the 
        ## corresponding records 
        df_seq = pd.read_csv("data/human.%s.upstream.score.and.filter.q.0.50.comb.txt"%chain, sep="\t", index_col=0)
        df_seq = df_seq[df_seq.scoreTag]  # filter to obtain only qualified records
        
        allele_seq_tuple_lst = zip(df_seq.index.tolist(), df_seq.Sequence.tolist())
        ref_seq_lst = [ "%s_%s"%(seq, allele) for allele, seq in allele_seq_tuple_lst ]
        
        # read the file containing both discovered sequences and the sample information
        df_comb = pd.read_csv("data/human_%s_upstream_combined_summary.txt"%chain, sep="\t")
        df_comb = df_comb[df_comb.seq_ref.isin(ref_seq_lst)]  # filter out the unqualified records
         
        # extract only the useful columns and rename them
        df_comb = df_comb[["SampleId", "Allele", "leader"]].rename(columns={"SampleId":"sample", "Allele":"allele"})
        
        # insert a column to record gene information
        df_comb["gene"] = df_comb.allele.apply(lambda x:x.split("*")[0])  
        
        # transform into a dictionary with keys being gene and value being unique leader list
        gene_leader_dict = df_comb[["sample", "gene", "leader"]].groupby("gene").leader.unique().to_dict()
        tuple_lst = sorted(gene_leader_dict.items(), key=lambda x:int(x[0][4:].split("-")[0].split("/")[0].split("D")[0]))
        gene_lst = [ x[0] for x in tuple_lst ]
        
        gene_a_lst = []
        gene_b_lst = []
        overlap_n_lst = []
        for i, (gene_a, leader_lst_a) in enumerate(tuple_lst):
            for j, (gene_b, leader_lst_b) in enumerate(tuple_lst):
                if i > j:
                    gene_a_lst.append(gene_a)
                    gene_b_lst.append(gene_b)
                    gene_set_a = set(leader_lst_a)
                    gene_set_b = set(leader_lst_b)
                    overlap_n_lst.append(len(gene_set_a & gene_set_b))
        df_gene_overlap = pd.DataFrame({"gene_a":gene_a_lst, "gene_b":gene_b_lst, "overlap":overlap_n_lst})
        df_mat = df_gene_overlap.pivot_table(index="gene_a", columns="gene_b", values="overlap").loc[gene_lst[1:], gene_lst[:-1]]
        
        # visualization
        plt.figure(figsize=(15,15))
        sns.heatmap(df_mat, cmap="Reds", vmin=0, vmax=3, square=True)
        plt.savefig("human_%s_gene_leader_overlap_heatmap.jpg"%chain, dpi=300)
        
        
if __name__ == "__main__":
    
    ## Run main function
    main()
