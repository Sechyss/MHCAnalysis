import numpy as np
import pandas as pd
import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from Bio.Seq import Seq

# %%
HLa_agregated = pd.read_table('/Users/u2176312/Downloads/hla_ligand_atlas/HLA_aggregated.tsv', sep='\t')

with open('/Users/u2176312/Downloads/hla_ligand_atlas/sequences_benign_human.fasta', 'a') as f:
    for index, row in HLa_agregated.iterrows():
        sequence = Seq(row['peptide_sequence'])
        fastaheader = str(row['peptide_sequence_id'])
        seq_record = SeqRecord(sequence, fastaheader, description='')
        SeqIO.write(seq_record, f, 'fasta')

# %% Filtering those hits with low percentage of identity

blastp_table_1 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                               'CSP/Humanrecogn/Malaria_Human_blastp.tsv', sep='\t', header=None)
filtered = blastp_table_1[blastp_table_1[2] >= 70]
list_of_matching_malaria_blast = set(filtered[0].tolist())
list_of_matching_human_blast = set(filtered[1].tolist())

blastp_table_2 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                               'CSP/Humanrecogn/Human_Malaria_blastp.tsv', sep='\t', header=None)
filtered_2 = blastp_table_2[blastp_table_2[2] >= 70]
list2_of_matching_human_blast = set(filtered_2[0].tolist())
list2_of_matching_malaria_blast = set(filtered_2[1].tolist())

matching_human_blast = list_of_matching_human_blast.intersection(list2_of_matching_human_blast)
matching_malaria_blast = list_of_matching_malaria_blast.intersection(list2_of_matching_malaria_blast)

# %% Comparison between the benign sequences and the NCBI peptides

blastp_table_4 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                               'CSP/Humanrecogn/Human_Malaria_NCBI_blastp.tsv', sep='\t', header=None)
filtered_4 = blastp_table_4[blastp_table_4[2] >= 70]
list4_of_matching_human_blast = set(filtered_4[0].tolist())
list4_of_matching_malaria_blast = set(filtered_4[1].tolist())
columns = ['query id', 'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens',
           'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']
filtered_4.columns = columns

# %%

tempfile = open('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/NCBI_CSP_peptides_11kmer.pickle', 'rb')
dictionary = pickle.load(tempfile)

Table = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                      'NCBI_CSP/NCBI_Pf3D7_blastp.tsv', sep='\t', header=None)
Table.columns = columns

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                        'CSP/NCBI_CSP/NCBI_CSP_peptides_11kmer_filtered.fasta', 'fasta')
dictionary_index = {}

for seq_record in fastafile:
    sequence = str(seq_record.seq)
    if sequence in dictionary.keys():
        dictionary_index.update({seq_record.id: dictionary[sequence]})

kmer_list = {}
dictionary_lengths = {}
for item in dictionary_index.keys():
    list_ncbi_sequences = dictionary_index[item]
    empty_list = []
    for sequence in list_ncbi_sequences:
        parts = sequence.rsplit('_')
        empty_list.append('_'.join(parts[:3]))

    kmer_list.update({item: list(set(empty_list))})
    dictionary_lengths.update({item: len(list_ncbi_sequences)})

filtered_4['Kmer position'] = np.nan
filtered_4['Number of NCBI sequences'] = np.nan

Table['Kmer position'] = np.nan
Table['Number of NCBI sequences'] = np.nan

for index, row in filtered_4.iterrows():
    key_to_dict = row['subject id']
    filtered_4['Kmer position'][index] = kmer_list[key_to_dict]
    filtered_4['Number of NCBI sequences'][index] = dictionary_lengths[key_to_dict]

for index, row in Table.iterrows():
    key_to_dict = row['subject id']
    Table['Kmer position'][index] = kmer_list[key_to_dict]
    Table['Number of NCBI sequences'][index] = dictionary_lengths[key_to_dict]
# filtered_4.to_csv('/Users/u2176312/OneDrive - University of Warwick/'
#                  'CSP/Humanrecogn/Blastp_NCBI_benign_location.csv')

Table.to_csv('/Users/u2176312/OneDrive - University of Warwick/'
             'CSP/Netchop/Blastp_NCBI_Pf3D7.csv')
