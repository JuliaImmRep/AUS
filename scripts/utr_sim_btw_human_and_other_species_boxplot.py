import pandas as pd
import csv, sys
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import pairwise2
from numpy import mean

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

def main():
    
    for i, species in enumerate(["rhesus", "fascicularis", "mouse", "rat"]):
    
        # import the discovered sequences
        df_human = pd.read_csv("data/human.igh.upstream.score.and.filter.q.0.50.comb.txt", sep="\t")
        df_other = pd.read_csv("data/%s.igh.upstream.score.and.filter.q.0.50.comb.txt"%species, sep="\t")

        # import the pairwise identity matrix to determine the 
        # core genes of other species
        df_mat = pd.read_csv("human.%s.ighv.pairwise.iden.txt"%species, sep="\t", index_col=0)
        mapDict = {}  # with key being the gene id of human and the value being a list of gene id of other species
        for alleleO, alleleH in df_mat.apply(lambda x:x.sort_values().index[-1], axis=1).to_dict().items():
            geneO = alleleO.split("*")[0]
            geneH = alleleH.split("*")[0]
            try:
                mapDict[geneH].append(geneO)
            except:
                mapDict.update({geneH: [geneO]})

        # obtain the sequence list for each gene for human
        df_human = df_human[df_human.scoreTag]
        df_human["gene"] = df_human.Allele.str.split("*", expand=True).iloc[:,0]
        df_human["utr"] = df_human.apply(lambda x:x.Sequence[:x.leaderStart], axis=1)
        df_human = df_human.groupby("gene").utr.unique()
        geneLeaderDictH = dict([(gene, list(seqLst)) for gene, seqLst in df_human.to_dict().items()])

        # obtain the sequence list for each gene for other species
        df_other = df_other[df_other.scoreTag]
        df_other["gene"] = df_other.Allele.str.split("*", expand=True).iloc[:,0]
        df_other["utr"] = df_other.apply(lambda x:x.Sequence[:x.leaderStart], axis=1)
        df_other = df_other[df_other.utr != ""]
        df_other = df_other.groupby("gene").utr.unique()

        geneLeaderDictO = df_other.to_dict()
        geneLeaderDictO = dict([(gene, list(seqLst)) for gene, seqLst in df_other.to_dict().items()])

        # obtain the core gene list 
        df_core = pd.read_csv("core_V_gene.txt", header=None)
        coreList = df_core[0].tolist()

        # calculate the leader similarity between core and noncore gene
        coreLst = []
        idenList = []
        for geneH in coreList:
            if geneH in geneLeaderDictH.keys():
                seqLstO = []
                flag = 0
                if geneH in mapDict.keys():
                    for geneO in mapDict[geneH]:
                        if geneO in geneLeaderDictO.keys():
                            if flag == 0:
                                coreLst.append(geneH)
                            seqLstO = seqLstO + geneLeaderDictO[geneO]
                            flag += 1
                if seqLstO != []:
                    seqLstH = geneLeaderDictH[geneH]
                    idenLst = []
                    for seqH in seqLstH:
                        for seqO in seqLstO:
                            #seq1, seq2 = pairwise2.align.globalxx(seqH, seqO)[0][:2]
                            seq1, seq2 = pairwise2.align.globalms(seqH, seqO, 2, -1, -5, -1)[0][:2]
                            idenLst.append(identity(seq1, seq2))
                    idenList.append(mean(idenLst))
        df_core = pd.DataFrame({"gene":coreLst, "identity":idenList})
        df_core["type"] = "core"
        df_core["species"] = species

        noncoreLst = list(set(geneLeaderDictH.keys()) - set(coreList))
        noncoreList = []
        idenList = []
        for geneH in noncoreLst:
            if geneH in geneLeaderDictH.keys():
                seqLstO = []
                flag = 0
                if geneH in mapDict.keys():
                    for geneO in mapDict[geneH]:
                        if geneO in geneLeaderDictO.keys():
                            if flag == 0:
                                noncoreList.append(geneH)
                            seqLstO = seqLstO + geneLeaderDictO[geneO]
                            flag += 1
                if seqLstO != []:
                    seqLstH = geneLeaderDictH[geneH]
                    idenLst = []
                    for seqH in seqLstH:
                        for seqO in seqLstO:
                            #seq1, seq2 = pairwise2.align.globalxx(seqH, seqO)[0][:2]
                            seq1, seq2 = pairwise2.align.globalms(seqH, seqO, 2, -1, -5, -1)[0][:2]
                            idenLst.append(identity(seq1, seq2))
                    idenList.append(mean(idenLst))
        df_noncore = pd.DataFrame({"gene":noncoreList, "identity":idenList})
        df_noncore["type"] = "noncore"
        df_noncore["species"] = species

        df = pd.concat([df_core, df_noncore])
        
        if i != 0:
            dF = pd.concat([dF, df])
        else:
            dF = df
    dF.to_csv("human_other_species_utr_similarity_core_and_noncore_igh.csv")

    df_data = pd.read_csv("human_other_species_utr_similarity_core_and_noncore_igh.csv", index_col=0)

    palette = ["limegreen", "tomato"]

    plt.figure(figsize=(4,3))
    ax = sns.boxplot(x="species", y="identity", hue="type", 
                     data=df_data, saturation=1,
                    palette=palette, hue_order = ["core", "noncore"],
                        order=["rhesus", "fascicularis", "mouse", "rat"])
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.set_ylim(0.25, 1.05)
    #ax.set_xticklabels("")
    #ax.set_yticklabels("")
    #ax.set_xlabel("")
    #ax.set_ylabel("")
    #ax.get_legend().remove()
    plt.savefig("utr_similarity_cmp_bwt_core_noncore_for_diff_species_igh.jpg", dpi=300, bbox_inches="tight") 

if __name__ == "__main__":
    main()
