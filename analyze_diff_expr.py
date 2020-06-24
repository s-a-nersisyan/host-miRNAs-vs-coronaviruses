import pandas as pd
import numpy as np

import random
from scipy.stats import hypergeom


#################
# Intersect DE miRNAs in two datasets and calculate p-value on intersection cardinality
#################

df1 = pd.read_csv("input_data/GSE36971/diff_expressed.csv", index_col=0)
df2 = pd.read_csv("input_data/GSE90624/diff_expressed.csv", index_col=0)

intersection = set.intersection(set(df1.index.tolist()), set(df2.index.tolist()))

df = df1.loc[intersection].join(df2.loc[intersection], lsuffix="1", rsuffix="2")
df = df[["baseMean1", "log2FoldChange1", "padj1", "baseMean2", "log2FoldChange2", "padj2"]]
# Remove miRNAs with different fold change sign
df = df.loc[df["log2FoldChange1"] / df["log2FoldChange2"] > 0]
common_miRNAs = df.index.to_list()
df = df.rename(columns={
    "baseMean1": "Control mean, GSE36971", 
    "log2FoldChange1": "Log2(fold change), GSE36971",
    "padj1": "Adjusted p-value, GSE36971",
    "baseMean2": "Control mean, GSE90624",
    "log2FoldChange2": "Log2(fold change), GSE90624",
    "padj2": "Adjusted p-value, GSE90624"
})
df.index.name = "pre-miRNA.miRNA"
df.to_csv("tables/mouse_lung_DE_miRNAs.tsv", sep="\t")

# Compute p-value using Monte-Carlo simulation
box = list(range(2302))  # There are 2303 miRNAs in total
N, counter = 10**2, 0  # FIXME
for i in range(N):
    set1 = set(random.sample(box, len(df1)))
    set2 = set(random.sample(box, len(df2)))
    if len(set.intersection(set1, set2)) >= len(df):
        counter += 1
print("p-value on two sets intersection: {}".format(counter / N))

#################
# Now analyze DE mRNA targets of miRNAs found
#################

df = pd.read_csv("input_data/GSE52405/diff_expressed.csv", index_col=0)
df = df.loc[np.abs(df["log2FoldChange"]) >= np.log2(1.5)]
# Add gene symbols
df_ens = pd.read_csv("input_data/Ensembl_mouse.tsv", index_col=0, sep="\t")
df = df.join(df_ens)

M = 19839  # Total number of genes
n = len(df.loc[df["log2FoldChange"] <= 0])  # Number of down-regulated genes

for miRNA in common_miRNAs:
    miRNA = miRNA.split(".")[1]
    
    df1 = pd.read_csv("input_data/miRDB/{}.tsv".format(miRNA), sep="\t", header=None)
    df2 = pd.read_csv("input_data/TargetScan/{}.tsv".format(miRNA), sep="\t")

    df1 = df1.loc[df1[1] >= 75]  # Target score >= 75
    df2 = df2.iloc[0:len(df2) // 10]  # Top-10% of TargetScan predictions (weighted context++ score)
    
    # Intersect TargetScan and miRDB
    genes1 = {g.strip() for g in df1[3].tolist()}
    genes2 = set(df2["Target gene"].tolist())
    common_TargetScan_miRDB = set.intersection(genes1, genes2)

    df_miRTarBase = pd.read_excel("input_data/miRTarBase8_mouse.xls")
    miRTarBase = {g.capitalize() for g in df_miRTarBase.loc[df_miRTarBase["miRNA"] == miRNA]["Target Gene"].to_list()}

    for gene_list, source in [(common_TargetScan_miRDB, "TargetScan + miRDB"), (miRTarBase, "miRTarBase")]:
        diff_expressed = [g for g in gene_list if g in df["Gene name"].to_list()]
        
        df_de = df.loc[df["Gene name"].isin(diff_expressed)]
        df_de = df_de.loc[df_de["log2FoldChange"] < 0]  # Only down-regulated
        diff_expressed = df_de["Gene name"].to_list()
        
        N = len(gene_list)
        rv = hypergeom(M, n, N)
        pval = sum([rv.pmf(k) for k in range(len(diff_expressed), N + 1)])
        print("miRNA: {}, source: {}. {} out of {} genes are down-regulated, p = {}".format(
            miRNA, source, len(diff_expressed), N, pval
        ))
