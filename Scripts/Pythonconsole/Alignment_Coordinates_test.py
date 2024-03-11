import os

import pandas as pd
import pickle
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

os.chdir('/Users/u2176312/OneDrive - University of Warwick/Otherproteins/Multiple_Alignment/')

alignment = AlignIO.read('Liverstageantigen1_Nterminal.muscle.faa', format="fasta")
counter = 1
coordinates = []
for record in alignment:
    if record.id == 'XP_024329106.1':
        sequence = record.seq
        for aminoacid in sequence:
            if aminoacid != '-':
                coordinates.append(counter)
                counter += 1
            else:
                coordinates.append(counter)
    else:
        continue

df = pd.DataFrame(index=[records.id for records in alignment], columns=coordinates)

print('Collecting coordinates for ' + 'Liverstageantigen1_Nterminal')
for i, col in tqdm(enumerate(alignment)):
    aligment_id = col.id
    all_seq = col.seq
    for index, aminoacid in enumerate(all_seq):
        df.loc[aligment_id, coordinates[index]] = aminoacid

Gene_depth = {}
print('Estimating the gene depth for ' + 'Liverstageantigen1_Nterminal')
ending_pos = 0
for position in range(len(df.columns)):
    if ending_pos != df.columns[-1]:
        starting_pos = df.columns[position]
        ending_pos = df.columns[position + 10]
        filtered_df = df[df.columns[position: position + 10]]
        count_rows_with_hyphen = filtered_df.eq('-').any(axis=1).sum()
        depth = (len(df.index) - count_rows_with_hyphen) / len(df.index)
        Gene_depth.update({starting_pos: depth})

Sequences_dictionary = {}

print('Writing Kmer sequences for ' + 'Liverstageantigen1_Nterminal')
with open(str('Liverstageantigen1_Nterminal') + '_kmers.fasta', 'a') as f:

    for index, row in tqdm(df.iterrows()):
        filtered_row = row[row != '-']
        all_coordinates = filtered_row.index.tolist()
        full_sequence = "".join(filtered_row.values.tolist())

        ending_pos = 0
        for position in range(len(all_coordinates)):
            if ending_pos != all_coordinates[-1]:
                starting_pos = all_coordinates[position]
                ending_pos = all_coordinates[position + 10]

                kmer = "".join(filtered_row.values[position: position + 11])
                id_kmer = "".join(['Liverstageantigen1_Nterminal', '_', str(starting_pos), '_', str(ending_pos), '_', index])
                if kmer not in Sequences_dictionary.keys():
                    Sequences_dictionary.update({kmer: [id_kmer]})
                else:
                    Sequences_dictionary[kmer].append(id_kmer)
                seq_record = SeqRecord(seq=Seq(kmer), id=id_kmer, description='')
                SeqIO.write(seq_record, f, 'fasta')
            else:
                continue


with open(str('Liverstageantigen1_Nterminal') + '.pickle', 'wb') as f2:
    pickle.dump(Sequences_dictionary, f2)

with open(str('Liverstageantigen1_Nterminal') + '_kmers_filtered.fasta', 'a') as f3:
    Sequence_counter = 0

    for key in tqdm(Sequences_dictionary.keys()):
        peptide = key
        sequence_id = 'Sequence_' + str(Sequence_counter)
        seq_record = SeqRecord(Seq(peptide), id=sequence_id, description='')
        SeqIO.write(seq_record, f3, 'fasta')
        Sequence_counter += 1


