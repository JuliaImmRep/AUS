'''
    This script is for spliting the gene into
    core gene and noncore gene.
'''
import pandas as pd
import csv, sys
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def main():
    # import the pairwise identity matrix to determine the 
    # core genes of other species
    for i, species in enumerate(["rhesus", "fascicularis", "mouse", "rat"]):
        df = pd.read_csv("human.%s.ighv.pairwise.iden.txt"%species, sep="\t", index_col=0)
        # determine the best aligned human gene
        map_series = df.apply(lambda x:x.sort_values().index[-1].split("*")[0], axis=1)
        map_series.name = "geneH"
        map_df = map_series.to_frame()
        df_mat = df.T.max().reset_index().rename(columns={"index":"allele", 0:"identity"})
        df_mat["species"] = species
        df_mat.index = df_mat.allele
        df_mat = df_mat.merge(map_df, right_index=True, left_index=True)
        if i != 0:
            df_comb = pd.concat([df_comb, df_mat])
        else:
            df_comb = df_mat
    
    df_comb= df_comb.groupby(["species", "geneH"]).identity.mean().reset_index()
    
    # annotate core gene obtain the core gene list 
    df_core = pd.read_csv("core_V_gene.txt", header=None)
    coreList = df_core[0].tolist()
    df_comb["type"] = df_comb.geneH.apply(lambda x:"core" if x in coreList else "noncore")
        
    # visualization the similarity
    plt.figure(figsize=(4,3))
    palette = ["limegreen", "tomato"]
    ax = sns.boxplot(x="species", y="identity", hue="type", 
                     palette=palette, saturation=1.0, 
                     data=df_comb, hue_order = ["core", "noncore"],
                    order=["rhesus", "fascicularis", "mouse", "rat"])
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.set_ylim(0.25, 1.05)
    ax.set_xticklabels("")
    ax.set_yticklabels("")
    ax.set_xlabel("")
    ax.set_ylabel("")
    #ax.set_yticks([0.30, 0.4, 0.50, 0.6, 0.70, 0.8, 0.90])
    ax.get_legend().remove()
    plt.savefig("v_gene_similarity_cmp_bwt_diff_species_igh.jpg", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
