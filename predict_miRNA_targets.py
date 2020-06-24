import os

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt


#################
# Intersect TargetScan and miRDB predictions and take miRNAs having miRDB Target score >= 75
#################

# Get list of all viruses
viruses = sorted(os.listdir("input_data/miRDB"))

# Intersect TargetScan and miRDB predictions
predicted_miRNA_dfs = []
num_summary = []  # Number of miRNAs predicted by each of tools
for virus in viruses:
    df_TargetScan = pd.read_csv("input_data/TargetScan/{}".format(virus), sep="\t")
    df_miRDB = pd.read_csv("input_data/miRDB/{}".format(virus), sep="\t")
    target_scores = {miRNA: score for miRNA, score in zip(df_miRDB["miRNA Name"], df_miRDB["Target Score"])}
    
    result_list = []
    TargetScan_miRNAs = set()  # This is used to count number of distinct miRNAs predicted by TargetScan
    for ind, row in df_TargetScan.iterrows():
        for name in row["miRNA_family_ID"].split("/"):
            # Add miR- prefix if name starts with digit
            new_name = "miR-" + name if name[0].isdigit() else name
            # Add hsa- prefix to match with miRDB
            new_name = "hsa-" + new_name

            TargetScan_miRNAs.add(new_name)
            
            if new_name not in target_scores:  # Intersect with miRDB
                continue

            target_score = target_scores[new_name]
            result_list.append([new_name, virus.rstrip(".txt"), row["MSA_start"], row["MSA_end"], row["Site_type"], target_score])
    
    # Resuling df
    df = pd.DataFrame(result_list, columns = ["miRNA", "Virus", "Start", "End", "Site type", "Target Score"])
    df.to_csv("TargetScan_miRDB_intersection/{}.tsv".format(virus.rstrip(".txt")), sep="\t", index=None)
    
    num_TargetScan = len(TargetScan_miRNAs)
    num_miRDB = len(df_miRDB)
    num_common = len(np.unique(df["miRNA"]))
    num_summary.append([virus.rstrip(".txt"), num_miRDB, num_TargetScan, num_common])

    # Filter by Target score >= 75
    best_miRNAs = df.loc[df["Target Score"] >= 75, ["miRNA", "Virus", "Target Score"]].drop_duplicates()
    predicted_miRNA_dfs.append(best_miRNAs)

# Save results
df = pd.concat(predicted_miRNA_dfs, axis=0)
df = df.pivot(index="miRNA", columns="Virus", values="Target Score")
df = df.dropna()  # Select miRNAs having target score >= 75 for each of viruses
df.to_csv("tables/miRNAs_vs_viruses_best_predictions.tsv", sep="\t")

df_num = pd.DataFrame(num_summary, columns=["Virus", "miRDB", "TargetScan", "Common"])
df_num.to_csv("tables/TargetScan_miRDB_number_of_predictions.tsv", sep="\t", index=None)

#################
# Plot heatmap based on miRNA-virus target scores
#################

df = pd.read_csv("tables/miRNAs_vs_viruses_best_predictions.tsv", sep="\t", index_col=0)
# List of high confidence miRNAs (according to miRBase) was manually intersected with the df
high_conf_miRNAs = [
    "hsa-miR-16-5p",
    "hsa-miR-195-5p",
    "hsa-miR-21-3p",
    "hsa-miR-3065-5p",
    "hsa-miR-421",
    "hsa-miR-424-5p"
]

df = df.loc[high_conf_miRNAs]
df = df.loc[df.mean(axis=1).sort_values(ascending=False).index]  # Sort by mean target score

p = sns.clustermap(
        df, row_cluster=False,
        cmap=sns.color_palette("Blues"), linecolor="lightgrey", linewidth=0.1,
        cbar_kws={"label": "Target score"}, dendrogram_ratio=(.05, .2)
)
plt.setp(p.ax_heatmap.get_xticklabels(), rotation=45)
p.ax_heatmap.set_xlabel("")
p.ax_heatmap.set_ylabel("")
plt.tight_layout()
plt.savefig("figures/target_scores.pdf")
plt.close()

#################
# Plot expression of miRNAs in human lungs (TCGA-LUAD)
#################

df = pd.read_csv("input_data/TCGA-LUAD_miRNA.tsv", sep="\t", index_col=0)
df = df.loc[high_conf_miRNAs].reset_index().rename(columns={"index": "miRNA"})
df = df.loc[df.mean(axis=1).sort_values(ascending=False).index]  # Sort by mean expression
df = pd.melt(df, id_vars="miRNA", var_name="Sample", value_name="Expression, log(RPM)")

cmap = sns.color_palette("Blues")
sns.boxplot(
        x="miRNA", y="Expression, log(RPM)", data=df, 
        color=sns.color_palette("Blues")[3], saturation=1
)
plt.xticks(rotation="45")
plt.xlabel("")
plt.tight_layout()
plt.savefig("figures/miRNA_TCGA-LUAD.pdf")
plt.close()
