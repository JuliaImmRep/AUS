'''
    This script is for calculating and visualizing
    5'UTR sequence similarity.
    
    Written by Xiujia Yang, on July 2020
'''

import pandas as pd
import seaborn as sns
from Bio import pairwise2
from Bio import SeqIO
from numpy import mean
import matplotlib.pyplot as plt

def identity(seq1, seq2):
    '''
        This function is for comparing the distance
        of two zipped sequences
    '''
    iden = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 == c2:
            iden += 1
    return iden/float(len(seq1))

for chain in ["igh", "igk", "igl"]:
    # import the discovered sequences
    df_seq = pd.read_csv("data/human.%s.upstream.score.and.filter.q.0.50.comb.txt"%chain, sep="\t", index_col=0)
    df_seq = df_seq[df_seq.scoreTag]  # filter to obtain only qualified records

    # obtain the unique leader sequence set for each gene
    df_seq["gene"] = df_seq.apply(lambda x: x.name.split("*")[0], axis=1)
    df_seq["utr"] = df_seq.apply(lambda x: x.Sequence[:x.leaderStart], axis=1)
    gene_leaders_dict = df_seq.groupby("gene")["utr"].unique().to_dict()

    # compared the leader similarity between gene
    iden_dict = {}
    for geneA, leaderSetA in gene_leaders_dict.items():
        leader_iden_lst = []
        gene_lst = []
        for geneB, leaderSetB in gene_leaders_dict.items():
            gene_lst.append(geneB)
            temp_iden_lst = []
            if geneA != geneB:
                for leaderA in leaderSetA:
                    for leaderB in leaderSetB:
                        seq1, seq2 = pairwise2.align.globalms(leaderA, leaderB, 2, -1, -5, -1)[0][:2]
                        temp_iden_lst.append(identity(seq1, seq2))
                leader_iden_lst.append(mean(temp_iden_lst))
            else:
                leader_iden_lst.append(0)
        iden_dict.update({geneA: leader_iden_lst})

    pd.DataFrame(iden_dict, index=gene_lst).to_csv("human.%s.utr.similarity.matrix.csv"%chain)


    # import the precalculated statistics
    df = pd.read_csv("human.%s.utr.similarity.matrix.csv"%chain, index_col=0)

    # sort index
    df = df.sort_index().T.sort_index().T
    gene_lst = df.index.tolist()
    fam_lst = [ int(x[4:].split("-")[0].split("/")[0].split("D")[0].split("S")[0]) for x in df.index ]
    gene_lst_sorted = [x[0] for x in sorted(zip(gene_lst, fam_lst), key=lambda x:x[1]) ]
    df = df.loc[gene_lst_sorted, gene_lst_sorted]

    # visualize
    gene_lst = df.index.tolist()
    fam_lst = [ x[4:].split("-")[0].split("/")[0].split("D")[0] for x in df.index ]
    if chain == "igh":
        color_map = {
            "1": (237, 125, 49),
            "2": (208, 206, 206),
            "3": (146, 208, 80),
            "4": (91, 155, 213),
            "5": (118, 113, 113),
            "6": (118, 113, 113),
            "7": (118, 113, 113),
            "8": (118, 113, 113),
            "9": (118, 113, 113),
            "10": (118, 113, 113)
        }
    else:
        color_map = {
        "1": (237, 125, 49),
        "2": (208, 206, 206),
        "3": (146, 208, 80),
        "4": (118, 113, 113),
        "5": (118, 113, 113),
        "6": (118, 113, 113),
        "7": (118, 113, 113),
        "8": (118, 113, 113),
        "9": (118, 113, 113),
        "10": (118, 113, 113)
        }
    color_maps = {}
    for k, v in color_map.items():
        v1 = v[0]/float(256)
        v2 = v[1]/float(256)
        v3 = v[2]/float(256)
        color_maps.update({k: (v1, v2, v3)})
    color_lst = [ color_maps[fam] for fam in fam_lst ]
    colors = pd.Series(color_lst, index=gene_lst)
    g = sns.clustermap(data=df, cmap="Reds", vmin=0, vmax=1,
                       row_cluster=False, col_cluster=False,
                   row_colors=colors, col_colors=colors)
    ax = g.ax_heatmap
    ax.set_xticks([])
    ax.set_yticks([])
    g.cax.remove()  # remove colorbar
    plt.savefig("human.%s.utr.similarity.matrix.jpg"%chain, dpi=50, bbox_inches="tight")
