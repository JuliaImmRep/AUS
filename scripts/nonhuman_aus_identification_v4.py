'''
    This script is for detection of V gene upstream 
    (utr+leader) sequences. 
    
    Written by Xiujia Yang, on July 7, 2020
'''

import csv, os, glob, sys
from Bio import SeqIO
from collections import Counter
from numpy import mean
import re
from scipy.stats import multinomial
from math import log
from scipy.stats import multinomial
from math import log

def genotyping(c1, c2):
    '''
        This function accepts the numbers of the
        top two alleles for a typical gene and 
        return the inferred genotype.
    '''
    # Define priors and make transformation to
    # account for the exceptions
    HH = [1, 0]; HH = [(x+0.01)/1.02 for x in HH]
    HD = [0.6,0.4]; HD = [(x+0.01)/1.02 for x in HD]
    # Make the input tuple and calculate the 
    # probability of different priors given 
    # the observed number of different alleles
    countTuple = (c1, c2)
    ph = multinomial.pmf(countTuple, n=sum(countTuple), p=HH)
    pd = multinomial.pmf(countTuple, n=sum(countTuple), p=HD)
    phLog = -1000 if ph == 0 else log(ph, 10)
    pdLog = -1000 if pd == 0 else log(pd, 10)
    # Determine the genotype according to the 
    # probability difference
    if phLog == -1000 and pdLog == -1000:
        return "Ambiguous"
    elif phLog == -1000:
        return "Heterozygote"
    elif pdLog == -1000:
        return "Homozygote"
    elif log(ph/pd, 10) >= 1:
        return "Homozygote"
    elif log(pd/ph, 10) >= 1:
        return "Heterozygote"
    else:
        return "Ambiguous"

def obtain_upstream(aligntag, seq, primer, allele, len_dict):
    '''
        This function is for obtaining the upstream sequence for 
        a given read. It accepts aligntags and target sequences,
        and output the upstream sequence and its associated identity.
    '''
    if "," not in aligntag:  # for merged sequences
        ss, se, _, qs = map(int, aligntag.split("|")[:4])
        if ss != 0:  # for sequence without upstream sequences (the alignment not starting from 0)
            return ["", -1, 0]
        else:  # for sequence with upstream sequences (the alignment not starting from 0)
            upstreamSeq, leader_start = extract_upstream(seq, primer, qs, allele, len_dict)
            identity = 1-len(list(re.findall(r'\d+', aligntag.split("|")[5])))/(float(se-ss)+len(re.findall(r'I\d+[ATCG]', aligntag)))
            return [upstreamSeq, identity, leader_start]
    else:
        aligntagLst = aligntag.split(",")
        seqLst = seq.split(",")
        if aligntagLst[0] == "":
            ss, se, _, qs = map(int, aligntagLst[1].split("|")[:4])
            if ss != 0:
                return ["", -1, 0]
            else:
                upstreamSeq, leader_start = extract_upstream(seqLst[1], primer, qs, allele, len_dict)
                identity = 1-len(list(re.findall(r'\d+', aligntagLst[1].split("|")[5])))/(float(se-ss)+len(re.findall(r'I\d+[ATCG]', aligntagLst[1])))
                return [upstreamSeq, identity, leader_start]
        else:
            ss1, se1, _, qs1 = map(int, aligntagLst[0].split("|")[:4])
            if ss1 != 0:
                if aligntagLst[1] != "":
                    ss2, se2, _, qs2 = map(int, aligntagLst[1].split("|")[:4])
                    if ss2 != 0:
                        return ["", -1, 0]
                    else:
                        upstreamSeq, leader_start = extract_upstream(seqLst[1], primer, qs2, allele, len_dict)
                        identity = 1-len(list(re.findall(r'\d+', aligntagLst[1].split("|")[5])))/(float(se2-ss2)+len(re.findall(r'I\d+[ATCG]', aligntagLst[1])))
                        return [upstreamSeq, identity, leader_start]
                else:
                    return ["", -1, 0]
            else:
                upstreamSeq, leader_start = extract_upstream(seqLst[0], primer, qs1, allele, len_dict)
                identity = 1-len(list(re.findall(r'\d+', aligntagLst[0].split("|")[5])))/(float(se1-ss1)+len(re.findall(r'I\d+[ATCG]', aligntagLst[0])))
                return [upstreamSeq, identity, leader_start]         
            

def extract_upstream(seq, primer, qs, allele, len_dict):
    '''
        This function is for extracting the upstream sequence.
        It accepts three parameters, namely the read sequence, the primer
        sequence and position in the read sequences that aligned against
        the first nucleotide of reference sequence, it output the upstream
        sequence.
    '''
    # define parameters for finding primers
    #slidingLen = 100  # the sliding length [0:primerLen], [1:primerLen] means a sliding length of 1
    #nMismatch = 1  # the maximum number of mismatches that are acceptable
    #nG = 3  # the least number of G that are consider the signal of RACE starting bases
    #nLen = 13  # the maximum distance to the end of primer that are accepted as the valid starting base hit
    
    upstream_seq = seq[:qs]
    ATG_loci = []
    for i in range(len(upstream_seq)-2):
        if upstream_seq[i:i+3] == "ATG":
            ATG_loci.append(i)
    if ATG_loci == []:
        return ["", ""]
    else:
        leader_start = sorted(zip(ATG_loci, 
                                  [abs(len_dict[allele]-(len(upstream_seq)-x)) for x in ATG_loci]), 
                              key=lambda x:x[1])[0][0]
        return (upstream_seq[leader_start:], 0)
        
def hmd4str(s1, s2):
    '''
        Calculating hamming distance
    '''
    if len(s1) != len(s2):
        raise ValueError("not equal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))    

def imgt_flag(allele, seq, dic):
    '''
        This function is for determining the imgt flag
        0 -> consistent with IMGT record
        1 -> disconsistent with IMGT record
        2 -> no available records for this allele in IMGT
    '''
    if allele not in dic.keys():
        return 2
    elif seq != dic[allele]:
        if seq == "":
            return 5
        elif seq in dic[allele]:
            return 3
        elif dic[allele] in seq:
            return 4
        else:
            return 1
    else:
        return 0

def substring(leader, utr):
    '''
        This is function is for extracting the 15 bp substring
        downstream and upstream ATG. The substring can be extracted
        only if the length for both UTR and leader greater than 15
        bp and 18 bp (5+13), respectively.
    '''
    if len(leader) >= 18 and len(utr) >= 15:
        return utr[-15:] + leader[:18]
    else:
        return ""
    
def main():
    ####################### File input module ######################
    ## Import the clones.txt file
    cloneHandle = csv.reader(open(cloneFl, "rU"), delimiter="\t")
    next(cloneHandle)  # skip the title line
    id_bestv_dict = {}  # looks like {"0": "IGHV9-1*02", "0": "IGHV11-2*02"}
    for line in cloneHandle:
        cloneId, cloneCount, bestVHit = line[0], float(line[1]), line[5].split("(")[0]
        if cloneCount >= cloneSizeFilter:
            id_bestv_dict.update({cloneId: bestVHit})
    bestv_freq_dict = Counter(id_bestv_dict.values())
    #qualified = [ bestv for bestv, freq in bestv_freq_dict.items() if freq > cloneNumFilter ]
    ## Import the alignments.txt file
    
    ## Store the utr sequence in a dictionary
    utr_seq_dict = {}
    for rec in SeqIO.parse(imgtUTRRef, "fasta"):
        allele = rec.id.split("|")[1]
        utr_seq_dict.update({allele: str(rec.seq).upper()})
    
    #### Determine the optimal leader sequence for each allele ######
    len_dict = {}
    seq_dict = {}  # exclude partial sequence for optimal leader length determination
    leader_seq_dict = {}  # include all known leader sequences including partial
    for rec in SeqIO.parse(imgtLeaderRef, "fasta"):
        allele = rec.id.split("|")[1]
        leader_seq_dict.update({allele: str(rec.seq).upper()})
        if "partial" not in rec.id:
            len_dict.update({allele: len(rec.seq)})
            seq_dict.update({allele: str(rec.seq).upper()})
    len_infer_dict = {}
    allele_infer_lst = list(set(bestv_freq_dict.keys()) - set(len_dict.keys()))
    for allele in allele_infer_lst:
        # Replace with the most prevalent length based on all allele corresponding to the same gene
        known_allele_lst = [ alt for alt in len_dict.keys() if alt.split("*")[0] == allele.split("*")[0] ]
        # Replace with the most prevalent length based on all allele corresponding to the genes of the same family
        if known_allele_lst == []:
            known_allele_lst = [ alt for alt in len_dict.keys() if alt.split("-")[0] == allele.split("-")[0] ]
        # Replace with the most prevalent length based on all allele
        if known_allele_lst == []:
            known_allele_lst = len_dict.keys()    
        optimal_len = sorted(Counter([len_dict[alt] for alt in known_allele_lst]).items(), key=lambda x:x[1], reverse=True)[0][0]
        len_infer_dict.update({allele: optimal_len})
    # add the infer length to the len_dict
    for allele, length in len_infer_dict.items():
        len_dict.update({allele: length})
    
    
    ##################### Genotype module #########################
    allele_lst = bestv_freq_dict.keys()
    gene_lst = list(set([bestv.split("*")[0] for bestv in bestv_freq_dict.keys()]))
    geno_dict = {}  # looks like {"IGHV1-1": "Haplotype", "IGHV1-2": "Diplotype"}
    for gene in gene_lst:
        # obtaining the list of tuples (allele, allele frequency) corresponding to this gene
        allele_freq_tuple_sorted = sorted([(allele, bestv_freq_dict[allele]) for allele in allele_lst if allele.split("*")[0] == gene], 
               key=lambda x:x[1], reverse=True)
        if len(allele_freq_tuple_sorted) != 1:
            c1, c2 = [ freq for allele, freq in allele_freq_tuple_sorted ][:2]
            temp_allele_lst = [ allele for allele, freq in allele_freq_tuple_sorted[:2] ]
            genotype = genotyping(c1, c2)
        else:
            c1, c2 = allele_freq_tuple_sorted[0][1], 0
            temp_allele_lst = [allele_freq_tuple_sorted[0][0]]
            genotype = genotyping(c1, c2)
        if genotype == "Homozygote":
            geno_dict.update({gene: [genotype, [temp_allele_lst[0]]]})
        elif genotype == "Heterozygote":
            geno_dict.update({gene: [genotype, temp_allele_lst]})
        else:
            geno_dict.update({gene: ["Heterozyote", [temp_allele_lst[0]]]})
    
    
    ###### Identification of upstream for each sequence ######
    clone_upstream_dict = {}  # looks like {"0":[("TTAGC", 0.87), ("TTAGC", 0.88)]}
    alignHandle = csv.reader(open(alignFl, "rU"), delimiter="\t")
    next(alignHandle)  # skip the title line
    for line in alignHandle:
        cloneId, bestVHit, bestVAlignment, targetSequences = line[1:]
        # discard sequences not assembled into clone
        if cloneId == "-1":
            continue
        # discard sequences with bestVHit different from its clone's
        if bestVHit != id_bestv_dict[cloneId]:
            continue
        # discard sequence assigned to an allele not in the genotype
        if bestVHit.split("*")[0] not in geno_dict.keys():
            continue
        # discard sequence assigned to an allele not in the genotype
        if bestVHit not in geno_dict[bestVHit.split("*")[0]][1]:
            continue
        if id_bestv_dict.has_key(cloneId):
            upstreamSeq, identity, leader_start = obtain_upstream(bestVAlignment, targetSequences, primer, bestVHit, len_dict)
            if upstreamSeq != "":
                try:
                    clone_upstream_dict[cloneId].append((upstreamSeq, identity, leader_start))
                except:
                    clone_upstream_dict.update({cloneId:[(upstreamSeq, identity, leader_start)]})
    
    ###### Identification of consensus for each clone ######
    clone_upstream_con_dict = {}  #clone_upstream_consensus_dict
    for cloneId, seq_iden_tuple in clone_upstream_dict.items():
        upstreamLst = [x[0] for x in seq_iden_tuple]
        upstream_leader_dict = dict(list(set([(x[0], x[2]) for x in seq_iden_tuple])))
        upstreamFreqDict = Counter(upstreamLst)
        #upstreamFreqTupleList = sorted(Counter(upstreamLst).items(), key=lambda x:x[1], reverse=Ture)
        #upstreamFreqIdenTupleList = []
        upstreamFreqIdenTupleList = [] 
        for upstream in upstreamFreqDict.keys():
            iden_lst = []
            for upstream_temp, iden, _ in seq_iden_tuple:
                if upstream == upstream_temp:
                    iden_lst.append(iden)
            upstreamFreqIdenTupleList.append((upstream, upstreamFreqDict[upstream], mean(iden_lst), upstream_leader_dict[upstream]))
        #### Consider only the top1 that comprising more than 50% of total sequences
        sortedUpstreamFreqIdenTupleList = sorted(upstreamFreqIdenTupleList, key=lambda x:(x[1],x[2]), reverse=True)
        if sortedUpstreamFreqIdenTupleList[0][1]/float(len(seq_iden_tuple)) > 0.5:
            try:
                clone_upstream_con_dict[id_bestv_dict[cloneId]].append(sortedUpstreamFreqIdenTupleList[0])
            except:
                clone_upstream_con_dict.update({id_bestv_dict[cloneId]: [sortedUpstreamFreqIdenTupleList[0]]})
        
    
    ###### Identification of Consensus sequence across clones ######
    f = open("%s.upstream.sequence.fasta"%sample, "wb")
    #fieldLst = ["Allele", "Group", "Flag", "Rank", "alleleCount", "novelCount", "Sequence", "seqFreq", "seqIdentity", "SampleId", "seq_ref", "IMGT"]
    fieldLst = ["Allele", "Group", "Flag", "Rank", "alleleCount", "novelCount", "Sequence", "leaderStart", 
                "seqFreq", "cloneFreq", "readFreq", "seqIdentity", "SampleId", "Isotype", "seq_ref", 
                "IMGT_leader", "IMGT_utr", "leader", "utr", "ATG_Substring"]
    f.write("\t".join(fieldLst)+"\n")
    for allele, seqFreqIdenTuple in clone_upstream_con_dict.items():
        #seqFreqIdenTuple = [ x for x in seqFreqIdenTuple if x[3] != 0 ]  # remove those clones without utr
        alleleCount = len(seqFreqIdenTuple)  # calculate the number of clones assigned with this allele
        if alleleCount >= 10:  # consider allele with a number of valid clones at least 10
            # determine the value for the fields listed by Yan Zhu
            geno = geno_dict[allele.split("*")[0]][0]
            novelCount = len([x[0] for x in seqFreqIdenTuple if x not in seq_dict.values()])
            # calculate the frequency of each unique upstream sequence
            # and the mean identity
            upstream_leader_dict = dict(list(set([(x[0], x[3]) for x in seqFreqIdenTuple])))
            upstreamLst = [x[0] for x in seqFreqIdenTuple]
            upstreamFreqDict = Counter(upstreamLst)
            upstreamFreqIdenTupleList = []
            for upstream in upstreamFreqDict.keys():
                iden_lst = []
                readNum = 0
                for upstream_temp, freq, iden, _ in seqFreqIdenTuple:
                    if upstream == upstream_temp:
                        iden_lst.append(iden)
                        readNum += freq
                upstreamFreqIdenTupleList.append((upstream, upstreamFreqDict[upstream], mean(iden_lst), readNum))
            sortedUpstreamFreqIdenTupleList = sorted(upstreamFreqIdenTupleList, key=lambda x:(x[1],x[2]), reverse=True)
            # generate variable for output
            sequence = sortedUpstreamFreqIdenTupleList[0][0]
            cloneFreq = sortedUpstreamFreqIdenTupleList[0][1]
            readFreq = sortedUpstreamFreqIdenTupleList[0][3]
            seqFreq = sortedUpstreamFreqIdenTupleList[0][1]/float(alleleCount)
            seqIdentity = sortedUpstreamFreqIdenTupleList[0][2]
            seq_ref = "%s_%s"%(sequence, allele)
            leader_start = upstream_leader_dict[sequence]
            leader = sequence[leader_start:]
            utr = sequence[:leader_start]
            IMGT_leader = imgt_flag(allele, leader, leader_seq_dict)
            IMGT_utr = imgt_flag(allele, utr, utr_seq_dict)
            ATG_Substring = substring(leader, utr)
            top1Line = [allele, geno, flag, 1, "NA", alleleCount, sequence, leader_start,
                        seqFreq, cloneFreq, readFreq, seqIdentity, sample, isotype, seq_ref,
                        IMGT_leader, IMGT_utr, leader, utr, ATG_Substring]
            top1Line = map(str, top1Line)
            if len(sortedUpstreamFreqIdenTupleList) >= 2:
                sequence = sortedUpstreamFreqIdenTupleList[1][0]
                cloneFreq = sortedUpstreamFreqIdenTupleList[1][1]
                readFreq = sortedUpstreamFreqIdenTupleList[1][3]
                seqFreq = sortedUpstreamFreqIdenTupleList[1][1]/float(alleleCount)
                seqIdentity = sortedUpstreamFreqIdenTupleList[1][2]
                seq_ref = "%s_%s"%(sequence, allele)
                leader_start = upstream_leader_dict[sequence]
                leader = sequence[leader_start:]
                utr = sequence[:leader_start]
                IMGT_leader = imgt_flag(allele, leader, leader_seq_dict)
                IMGT_utr = imgt_flag(allele, utr, utr_seq_dict)
                ATG_Substring = substring(leader, utr)
                top2Line = [allele, geno, flag, 2, "NA", alleleCount, sequence, leader_start,
                            seqFreq, cloneFreq, readFreq, seqIdentity, sample, isotype, seq_ref, 
                            IMGT_leader, IMGT_utr, leader, utr, ATG_Substring]
                top2Line = map(str, top2Line)
            if len(sortedUpstreamFreqIdenTupleList) == 1:
                f.write("\t".join(top1Line)+"\n")
            else:
                if upstreamFreqDict[sortedUpstreamFreqIdenTupleList[0][0]]/float(alleleCount) >= 0.875:
                    f.write("\t".join(top1Line)+"\n")
                else:
                    if geno == "Homozygote":
                        f.write("\t".join(top1Line)+"\n")
                        f.write("\t".join(top2Line)+"\n")
                    else:
                        f.write("\t".join(top1Line)+"\n")
    
if __name__ == "__main__":
    cloneFl = sys.argv[1]
    alignFl = sys.argv[2]
    primer = sys.argv[3]
    sample = sys.argv[4]
    imgtLeaderRef = sys.argv[5]
    imgtUTRRef = sys.argv[6]
    isotype = sys.argv[7]
    cloneSizeFilter = 1
    flag = "upstream"
    ## Run main function
    main()
