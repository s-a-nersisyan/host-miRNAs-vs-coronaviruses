import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.spatial.distance import squareform
from scipy.cluster import hierarchy

from itertools import combinations


def print_alignment(left, right, padding):
    print(" "*17 + " "*padding + "|"*(right - left))
    for virus in MAS:
        if left in local_coordinates[virus]:
            local_start = local_coordinates[virus][left] + 1
        else:
            local_start = None
        if right in local_coordinates[virus]:
            local_end = local_coordinates[virus][right]
        else:
            local_end = None
        
        region = annotations[virus].get(local_start)

        print("{:17s}".format(virus), end="")
        print(MAS[virus][left - padding:right + padding], end=" ")
        print(local_start, local_end, region)
    print(" "*17 + " "*padding + "|"*(right - left))


def calc_matches(left, right, padding, viruses):
    matches = 0
    for i in range(left - padding, right + padding):
        nucl = {MAS[virus][i] for virus in viruses}
        if len(nucl) == 1:
            matches += 1

    return matches


def calc_matching_proba(padding, matches, viruses):
    N = 10000
    hits = 0
    for i in range(N):
        a = np.random.randint(padding, len(MAS["SARS-CoV"]) - padding - 6)
        b = a + 6
        if calc_matches(a, b, padding, viruses) >= matches:
            hits += 1
    
    return hits / N

#################
# Parse Multiple Sequence Alignment (MSA) of coronavirus genomes
# and establish mapping between local and global coordinates
#################

f = open("input_data/MAS_coronaviruses.clustal")
f.readline()
f.readline()

MAS = {}
num_viruses = 7
while 1:
    line = f.readline()  # Empty line or EOF
    if len(line) == 0:  # EOF
        break
   
    for i in range(num_viruses):
        line = f.readline().strip()
        virus, seq = line.split(" "*6)
        MAS[virus] = MAS.get(virus, "") + seq
    
    f.readline()  # Line with meta information
 
MAS_coordinates = {}  # local coordinates -> MAS coordinates
for virus in MAS:
    MAS_coordinates[virus] = {}
    i = 0
    for j, s in enumerate(MAS[virus]):
        if s == "-":
            continue
        MAS_coordinates[virus][i] = j
        i += 1

# Reverse MAS_coordinates to get MAS coordinates -> local coordinates
local_coordinates = {virus: {v: k for k, v in MAS_coordinates[virus].items()} for virus in MAS_coordinates}

#################
# Parse fasta files with genome annotations
#################

annotations = {}
for virus in MAS:
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
# Load binding site predictions for hsa-miR-21-3p and map them onto alignment
#################

miRNA = "hsa-miR-21-3p"
binding_positions = {}  # Interval (start, end) -> list of viruses
for virus in MAS:
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
        start_aln, end_aln = MAS_coordinates[virus][start - 1], MAS_coordinates[virus][end]
        
        #dict_ = binding_positions.get((start_aln, end_aln), {})
        #dict_[virus] = (start, end, annotation[virus].get(start, []))
        binding_positions[(start_aln, end_aln)] = binding_positions.get((start_aln, end_aln), []) + [virus]


#################
# Draw virus similarity dendrogram based on common miRNA binding sites
#################

# Number of common binding positions for two viruses
df_common = pd.DataFrame(np.zeros((len(MAS), len(MAS))))
df_common.index = sorted(list(MAS.keys()))[::-1]
df_common.columns = sorted(list(MAS.keys()))[::-1]

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

padding = 50
for (a, b), viruses in sorted(binding_positions.items(), key=lambda x: len(x[1]), reverse=True):
    if len(viruses) == 1:
        continue

    print("Consensus between " + ", ".join(viruses) + ":\n")
    print_alignment(a, b, padding)
    print("\n" + "*"*17 + "\n")
