'''
    This script performs statistical analyses of the 
    length of 5'UTR in uORF group and non-uORF group.
    
    Written by Xiujia Yang, on July 2020.
'''

import pandas as pd
import os
from Bio import SeqIO
import seaborn as sns
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

upstream_lst = []
leader_lst = []
utr_lst = []
df_lst = []  # for transformed into dataframe

i = 0 
for species in ["human", "fascicularis", "rhesus", "mouse", "rat"]:
    for chain in ["igh", "igk", "igl"]:
        seq_fl = "data/%s.%s.upstream.score.and.filter.q.0.50.comb.txt"%(species, chain)
        if not os.path.exists(seq_fl):
            continue
        else:
            i += 1
        # read the discovered sequences
        df = pd.read_csv(seq_fl, sep="\t", index_col=0)
        # remove those not qualified records
        df = df[df.scoreTag]
        # remove utr-absent upstream sequences (for rhesus heavy chain)
        df = df[df.leaderStart!=0]
        
        # made columns for UTR sequences and tag to indicate
        # if it is uAUG-containing
        df["UTR"] = df.apply(lambda x:x["Sequence"][:x.leaderStart], axis=1)
        df["Length"] = df.apply(lambda x:len(x.UTR), axis=1)
        df["Group"] = df.apply(lambda x:"uAUG-containing" if "ATG" in x["UTR"] else "uAUG-absent", axis=1)
        df["Species-Chain"] = "%s-%s"%(species, chain)
        # combine the dataset from different species-chain
        df_temp = df.loc[:,["Species-Chain", "Group", "Length"]]
        if i != 1:
            dF = pd.concat([dF, df_temp])
        else:
            dF = df_temp
        # calculate the p-value 
        lst1 = df_temp[df_temp.Group=="uAUG-containing"]["Length"].tolist()
        lst2 = df_temp[df_temp.Group=="uAUG-absent"]["Length"].tolist()
        n1 = len(lst1)
        n2 = len(lst2)
        p_value = ttest_ind(lst1,lst2)[1]
        print ("%s %s %d,%d,%f"%(species, chain, n1, n2, p_value))
    
    # visualize the length distributions
    dF = dF[dF["Species-Chain"] != "rhesus-igh"]
    order = ["human-igh", "human-igk", "human-igl", "rhesus-igk",
                "fascicularis-igh", "mouse-igh", "mouse-igk", "rat-igh", "rat-igk"]
    hue_order = ["uAUG-absent", "uAUG-containing"]
    ax = sns.barplot(x="Species-Chain", y="Length", hue="Group", data=dF, 
                      linewidth=0.5, palette=["#0099CC", "#FF6666"], order=order, hue_order=hue_order)
    ax.get_legend().remove()  # remove the long leader legend
    ax.set_xlabel("")
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.xaxis.set_tick_params(width=2)
    ax.yaxis.set_tick_params(width=2)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.xticks(rotation=45, ha='right')
    plt.ylim(0, 88)
    plt.savefig("length_cmp_bwt_AUG_containing_absent.jpg", dpi=800, bbox_inches="tight")
