'''
    This script is for investigating the correlation 
    between the identity of variable genes of different species
    and the identity of their associated leader sequences.
'''


import pandas as pd
import seaborn as sns
from Bio import pairwise2
from numpy import mean
from sklearn import linear_model
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def leaderIden(geneO, geneH, dicO, dicH):
    '''
        This function is for calculating the identity
        between leader.
    '''
    seqLstO = dicO[geneO]
    seqLstH = dicH[geneH]
    idenLst = []
    for seqO in seqLstO:
        for seqH in seqLstH:
            seq1, seq2 = pairwise2.align.globalms(seqH, seqO, 2, -1, -5, -1)[0][:2]
            idenLst.append(identity(seq1, seq2))
    return mean(idenLst)

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

for species in ["rhesus", "fascicularis", "mouse", "rat"]:
    # import the discovered sequences
    df_human = pd.read_csv("data/human.igh.upstream.score.and.filter.q.0.50.comb.txt", sep="\t")
    df_other = pd.read_csv("data/%s.igh.upstream.score.and.filter.q.0.50.comb.txt"%species, sep="\t")
    # obtain the leader sequence identity 
    df_human = df_human[df_human.scoreTag]
    df_human["leader"] = df_human.apply(lambda x:x.Sequence[x.leaderStart:], axis=1)
    df_human = df_human.groupby("Allele").leader.unique()
    geneLeaderDictH = dict([(gene, list(seqLst)) for gene, seqLst in df_human.to_dict().items()])
    df_other = df_other[df_other.scoreTag]
    df_other["leader"] = df_other.apply(lambda x:x.Sequence[x.leaderStart:], axis=1)
    df_other = df_other.groupby("Allele").leader.unique()
    geneLeaderDictO = dict([(gene, list(seqLst)) for gene, seqLst in df_other.to_dict().items()])
    
    # import the pairwise identity matrix and determine the map relation
    df_mat = pd.read_csv("human.%s.ighv.pairwise.iden.txt"%species, sep="\t", index_col=0)
    
    df_data =  df_mat.stack().reset_index().rename(columns={"level_0":"geneO", "level_1":"geneH", 0:"vGeneIden"})
    
    df_data = df_data[df_data.geneO.isin(geneLeaderDictO.keys()) & df_data.geneH.isin(geneLeaderDictH.keys())]
    df_data["leaderIden"] = df_data.apply(lambda x:leaderIden(x.geneO, x.geneH, geneLeaderDictO, geneLeaderDictH), axis=1).tolist()
    df_data.reset_index().to_csv("human.%s.v.leader.identity.csv"%species)

fig, axes = plt.subplots(1,4, figsize=(11,3))
for i, species in enumerate(["rhesus", "fascicularis", "mouse", "rat"]):
    df_res = pd.read_csv("human.%s.v.leader.identity.csv"%species, index_col=0)
    df_res["leaderIden"] = df_res["leaderIden"] * 100
    df_res["vGeneIden"] = df_res["vGeneIden"] * 100

    # scatterplot visualization
    ax = sns.regplot(x="vGeneIden", y="leaderIden", data=df_res, 
                     ax=axes[i], ci=90, fit_reg=True, 
                     scatter_kws={"s": 15, "edgecolors":"white", "linewidth":0.3}, color="royalblue")
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)

    # calculate regression equation
    regr = linear_model.LinearRegression()
    regr.fit(df_res.loc[:,'leaderIden'].astype(float).values.reshape(-1,1), df_res.loc[:,"vGeneIden"].astype(float))
    [coef], intercept = regr.coef_, regr.intercept_
    r2 =regr.score(df_res.loc[:,'leaderIden'].astype(float).values.reshape(-1,1), df_res.loc[:,"vGeneIden"].astype(float))
    
    ax.set_xlim(45, 95)
    ax.set_ylim(30, 105)
    ax.set_xticklabels("")
    ax.set_yticklabels("")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout(pad=3)
    # print the model parameters and p-value
    print (coef, intercept, r2, pearsonr(df_res["vGeneIden"], df_res["leaderIden"]))
plt.savefig("human.other.species.v.leader.identity.corr.jpg", dpi=300, bbox_inches="tight")
