#! /home/alexmakh/miniconda3/bin/python3

import pandas as pd

virus_names = ['MERS-CoV', 'HCoV-HKU1',
               'HCoV-NL63', 'HCoV-OC43', 'HCoV-229E']

virus_ids = {
    'MERS-CoV': 'NC_019843.3',
    'HCoV-HKU1': 'NC_006577.2',
    'HCoV-229E': 'NC_002645.1',
    'HCoV-NL63': 'NC_005831.2',
    'HCoV-OC43': 'NC_006213.1',
}

#################
# Parse Multiple Sequence Alignment (MSA) of coronavirus genomes
# and establish mapping between local and global coordinates
#################


for n in virus_names:
    f = open("input_data/MSA_genomes/aligned_{}.fasta".format(n))

    MAS = {}
    genomes = f.read().split('>')[1:]
    for i in genomes:
        lines = i.split('\n')
        header = lines[0]
        virus_id = header.split(' |')[0]
        sequence = ''.join(lines[1:])
        MAS[virus_id] = sequence
    print(len(MAS))

    MAS_coordinates = {}  # local coordinates -> MAS coordinates

    virus = virus_ids[n]
    MAS_coordinates[virus] = {}
    i = 0
    for j, s in enumerate(MAS[virus]):
        if s == "-":
            continue
        MAS_coordinates[i] = j
        i += 1

#################
# Load binding site predictions for hsa-miR-21-3p and
# count mutations in genomes
#################
    print('\n', 'coronavirus', n, '\n')
    miRNA = "hsa-miR-21-3p"
    df = pd.read_csv("TargetScan_miRDB_intersection/{}.tsv".format(n),
                     sep="\t", index_col=0)
    binding_positions = {}  # Interval (start, end) -> list of viruses
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
        start_aln, end_aln = MAS_coordinates[start - 1], MAS_coordinates[end]

        key_m = MAS[virus][start_aln:end_aln + 1]
        count_mut = 0
        for na in MAS:
            if MAS[na][start_aln:end_aln + 1] == key_m:
                continue
            count_mut += 1
        print('binding postion {}-{}:  '.format(start, end), count_mut)
