import pandas as pd
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
blastp_table_3 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                               'CSP/Humanrecogn/Malaria_NCBI_Human_blastp.tsv', sep='\t', header=None)
filtered_3 = blastp_table_3[blastp_table_3[2] >= 70]
list3_of_matching_malaria_blast = set(filtered_3[0].tolist())
list3_of_matching_human_blast = set(filtered_3[1].tolist())

blastp_table_4 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                               'CSP/Humanrecogn/Human_Malaria_NCBI_blastp.tsv', sep='\t', header=None)
filtered_4 = blastp_table_4[blastp_table_4[2] >= 70]
list4_of_matching_human_blast = set(filtered_4[0].tolist())
list4_of_matching_malaria_blast = set(filtered_4[1].tolist())

matching_human_blast_2 = list3_of_matching_human_blast.intersection(list4_of_matching_human_blast)
matching_malaria_blast_2 = list3_of_matching_malaria_blast.intersection(list4_of_matching_malaria_blast)