import pandas as pd
import numpy as np


def mutation_rate(start, end, canonical):
    '''
    Total number of mutations divided by total number of genomes and sequence length
    '''
    seq_stack = [MSA[v][start:end] for v in MSA]
    total_mutations = 0
    for i in range(len(canonical)):
        total_mutations += len([1 for seq in seq_stack if seq[i] != canonical[i]])
    
    return total_mutations / len(seq_stack) / len(canonical)


#################
# Parse Multiple Sequence Alignment (MSA) of coronavirus genomes
# and establish mapping between local and global coordinates
#################

virus_ids = {
    "SARS-CoV-2": "NC_045512.2",
    'SARS-CoV': 'NC_004718.3',
    'MERS-CoV': 'NC_019843.3',
    'HCoV-229E': 'NC_002645.1',
    'HCoV-HKU1': 'NC_006577.2',
    'HCoV-NL63': 'NC_005831.2',
    'HCoV-OC43': 'NC_006213.1',
}

mutation_rates = {}
for virus in virus_ids:
    f = open("input_data/MSA/{}.fasta".format(virus))

    MSA = {}
    genomes = f.read().split('>')[1:]
    for g in genomes:
        lines = g.split('\n')
        header = lines[0]
        virus_id = header.split(' |')[0]
        sequence = ''.join(lines[1:])
        
        if "N" in sequence:
            continue

        MSA[virus_id] = sequence
    
    MSA_coordinates = {}  # local coordinates -> MSA coordinates
    i = 0
    for j, s in enumerate(MSA[virus_ids[virus]]):
        if s == "-":
            continue
        MSA_coordinates[i] = j
        i += 1
    

    #################
    # Load binding site predictions for hsa-miR-21-3p and
    # count mutations in putative regions
    #################

    miRNA = "hsa-miR-21-3p"
    canonical_seed = "GGTGTT"
    df = pd.read_csv("TargetScan_miRDB_intersection/{}.tsv".format(virus), sep="\t", index_col=0)
    
    mutation_rates[virus] = {}
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
        start_aln, end_aln = MSA_coordinates[start - 1], MSA_coordinates[end]
        
        mutation_rates[virus][(start, end)] = mutation_rate(start_aln, end_aln, canonical_seed)

# Summary on mutation rates
for virus in mutation_rates:
    print("{}: number of binding sites without mutations: {} out of {}".format(
        virus, len([1 for v in mutation_rates[virus].values() if v == 0]), len(mutation_rates[virus])
    ))
    print("{}: mean average mutation rate: {}\n".format(
        virus, np.mean([v for v in mutation_rates[virus].values()])
    ))

#################
# Mutation rates for shared binding sites
#################

for line in open("tables/binding_sites.tsv"):
    if line.startswith("Virus"):
        pass
    elif line.strip() == "":
        print()
    else:
        virus, start, end, seq = line.strip().split("\t")
        print(virus, start, end, mutation_rates[virus][(int(start), int(end))])
