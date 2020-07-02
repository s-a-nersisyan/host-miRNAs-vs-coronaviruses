import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy

from itertools import combinations


def print_alignment(left, right, viruses, padding):
    to_print = ["Virus", "Seed start", "Seed end", "Alignment"]
    print("\t".join(to_print))
    for virus in viruses:
        if left in local_coordinates[virus]:
            local_start = local_coordinates[virus][left] + 1
        else:
            local_start = None
        if right in local_coordinates[virus]:
            local_end = local_coordinates[virus][right]
        else:
            local_end = None
        
        region = annotations[virus].get(local_start)
        
        to_print = [
            virus, str(local_start), str(local_end), 
            MSA[virus][left - padding:left] + "|" + MSA[virus][left:right] + "|" + MSA[virus][right:right + padding]
        ]
        print("\t".join(to_print))

#################
# Parse Multiple Sequence Alignment (MSA) of coronavirus genomes
# and establish mapping between local and global coordinates
#################

MSA = {}
virus_ids = {
    "NC_002645.1": "HCoV-229E",
    "NC_006577.2": "HCoV-HKU1",
    "NC_005831.2": "HCoV-NL63",
    "NC_006213.1": "HCoV-OC43",
    "NC_019843.3": "MERS-CoV",
    "NC_045512.2": "SARS-CoV-2",
    "NC_004718.3": "SARS-CoV",
}
f = open("input_data/MSA/human_coronaviruses.fasta")
genomes = f.read().split(">")[1:]
for g in genomes:
    lines = g.split("\n")
    header = lines[0]
    sequence = "".join(lines[1:])
    virus_id = header.split(" |")[0]
    virus = virus_ids[virus_id]
    MSA[virus] = sequence

MSA_coordinates = {}  # local coordinates -> MSA coordinates
for virus in MSA:
    MSA_coordinates[virus] = {}
    i = 0
    for j, s in enumerate(MSA[virus]):
        if s == "-":
            continue
        MSA_coordinates[virus][i] = j
        i += 1

# Reverse MSA_coordinates to get MSA coordinates -> local coordinates
local_coordinates = {virus: {v: k for k, v in MSA_coordinates[virus].items()} for virus in MSA_coordinates}

#################
# Parse fasta files with genome annotations
#################

annotations = {}
for virus in MSA:
    annotations[virus] = {}
    f = open("input_data/genomes/{}/cds.fasta".format(virus))
    for line in f:
        if not line.startswith(">"):
            continue
        

        if "join" not in line:
            name = line.split(" |")[1].split(" [")[0]
            region = line.split(" |")[0].split(":")[1]
            start, end = region.split("..")
        else:
            # Polypeptides
            name = line.split(")|")[1].split(" [")[0]
            
            region1 = line.split(" |")[0].split(":")[1]
            start = region1.split("..")[0]
            
            region2 = line.split(")|")[0].split(":")[-1]
            end = region2.split("..")[1]

        for i in range(int(start), int(end) + 1):
            annotations[virus][i] = annotations[virus].get(i, []) + [name]

#################
# Create summary table on hsa-miR-21-3p binding positions
#################

miRNA = "hsa-miR-21-3p"
table = []
for virus in sorted(MSA)[::-1]:
    df = pd.read_csv("TargetScan_miRDB_intersection/{}.tsv".format(virus), sep="\t", index_col=0)
    site_types = np.unique(df["Site type"]).tolist()
    nums = [len(df.loc[df["Site type"] == st].loc[[miRNA]]) for st in site_types]
    table.append([virus] + [str(n) for n in nums] + [str(sum(nums))])

# Print latex table (Table 1)
print("& " + " & ".join(site_types + ["Total"]) + " \\\\\\hline")
for row in table:
    print(" & ".join(row) + " \\\\\\hline")


#################
# Load binding site predictions for hsa-miR-21-3p and map them onto alignment
#################

miRNA = "hsa-miR-21-3p"
binding_positions = {}  # Interval (start, end) -> list of viruses
for virus in MSA:
    df = pd.read_csv("TargetScan_miRDB_intersection/{}.tsv".format(virus), sep="\t", index_col=0)
    for i, row in df.loc[[miRNA]].iterrows():
        start, end = row["Start"], row["End"]
        # Get classical 6-mer seed binding coordinates
        if row["Site type"] == '7mer-m8':
            start += 1
        elif row["Site type"] == '7mer-1a':
            end -= 1
        elif row["Site type"] == "8mer-1a":
            start += 1
            end -= 1
        
        # Transform to alignment coordinates
        start_aln, end_aln = MSA_coordinates[virus][start - 1], MSA_coordinates[virus][end]
        
        binding_positions[(start_aln, end_aln)] = binding_positions.get((start_aln, end_aln), []) + [virus]


#################
# Draw virus similarity dendrogram based on common miRNA binding sites
#################

# Number of common binding positions for two viruses
df_common = pd.DataFrame(np.zeros((len(MSA), len(MSA))))
df_common.index = sorted(list(MSA.keys()))[::-1]
df_common.columns = sorted(list(MSA.keys()))[::-1]

total_positions = 0
for (a, b), viruses in binding_positions.items():
    if len(viruses) == 1:
        continue
    
    for v1, v2 in combinations(sorted(viruses), 2):
        # For matrix simmetricity 
        df_common.loc[v1, v2] += 1
        df_common.loc[v2, v1] += 1
    
    total_positions += 1

mpl.rcParams['axes.titlepad'] = 20
plt.figure(figsize=(6, 6))
sns.heatmap(
    df_common, annot=True, annot_kws={"color": "white"}, cbar=False,
    cmap=sns.color_palette("Blues")[2:], linecolor="lightgrey", linewidth=0.1
)
plt.xticks(rotation="45")
plt.tick_params(axis="y", right=True, left=False, labelright=True, labelleft=False, labelrotation=0)
plt.title("A", loc="left", fontdict={"fontsize": "xx-large", "fontweight": "bold"})
plt.tight_layout()
plt.savefig("figures/shared_binding_sites.pdf")
plt.close()

# Convert to distance matrix
for v1, v2 in combinations(sorted(df_common.columns), 2):
    df_common.loc[v1, v2] = 1 - df_common.loc[v1, v2] / total_positions
    df_common.loc[v2, v1] = 1 - df_common.loc[v2, v1] / total_positions

plt.figure(figsize=(6, 6))

hierarchy.dendrogram(
        hierarchy.linkage(squareform(df_common), 'average'),
        labels=df_common.columns.to_list(), 
        color_threshold=0, 
        above_threshold_color='black',
        leaf_font_size=10
)
plt.xticks(rotation="45")
plt.title("B", loc="left", fontdict={"fontsize": "xx-large", "fontweight": "bold"})
plt.tight_layout()
plt.savefig("figures/binding_site_dendrogram.pdf")

#################
# Print results
#################

padding = 20
for (a, b), viruses in sorted(binding_positions.items(), key=lambda x: len(x[1]), reverse=True):
    if len(viruses) == 1:
        continue

    print_alignment(a, b, viruses, padding)
    print()
