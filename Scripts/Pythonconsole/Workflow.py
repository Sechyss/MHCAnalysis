import pandas as pd
from Bio import SeqIO

columns = ['query id', 'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens',
           'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']
Table_blastp = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                             'NCBI_CSP/NCBI_Pf3D7_blastp.tsv', sep='\t', header=None)
Table_blastp.columns = columns
Table_blastp = Table_blastp[Table_blastp['% identity'] >= 80]  # Filter to 80% identity

Table_blastp['subject id'] = Table_blastp['subject id'].apply(lambda x: int(str(x).replace('Sequence_', '')) + 1)
Highpercentage = list(set(Table_blastp['subject id'].tolist()))
Table_mhc = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                          'NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
                          'NCBI_TopABC_all_lengths_NCBIseqs.txt', sep='\t')
Table_mhc = Table_mhc[Table_mhc['seq_num'].isin(Highpercentage)]  # Filter only those sequences with matches in humans

mhc_dict = {}
for index, row in Table_mhc.iterrows():
    if row['seq_num'] not in mhc_dict.keys():
        mhc_dict.update({row['seq_num']: [row['allele']]})
    else:
        mhc_dict[row['seq_num']].append(row['allele'])

for key in mhc_dict.keys():
    new_data = list(set(mhc_dict[key]))
    mhc_dict[key] = new_data

sequence = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
                       'NCBI_CSP_peptides_11kmer_filtered.fasta', 'fasta')

dict_seq_C_terminal = {}

for seq_record in sequence:
    cterminal = str(seq_record.seq)[-1]
    sequenceID = int(str(seq_record.id).replace('Sequence_', ''))+1
    dict_seq_C_terminal.update({sequenceID: cterminal})


NetchopTable = pd.read_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/Netchop/result_CSP_Pf3D7_netchop.csv')
NetchopTable = NetchopTable[NetchopTable['prediction_score'] >= 0.5]  # Filter cuts with < 0.5 hits

netchop_dict = {}
for index, row in NetchopTable.iterrows():
    netchop_dict.update({row['#']: row['amino_acid']})

Table_blastp_Pf3D7 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                                   'NCBI_CSP/NCBI_Pf3D7_blastp.tsv', sep='\t', header=None)
Table_blastp_Pf3D7.columns = columns
Table_blastp_Pf3D7 = Table_blastp_Pf3D7[Table_blastp_Pf3D7['% identity'] >= 80]
Table_blastp_Pf3D7['subject id'] = Table_blastp_Pf3D7['subject id'].apply(lambda x: int(str(x).replace('Sequence_', '')) + 1)
Table_blastp_Pf3D7 = Table_blastp_Pf3D7[Table_blastp_Pf3D7['subject id'].isin(Highpercentage)]


Table_blastp_Pf3D7['Absolute start'] = ''
Table_blastp_Pf3D7['Absolute end'] = ''
Table_blastp_Pf3D7['HLA recognition'] = ''
Table_blastp_Pf3D7['C terminal chop'] = ''
for index, row in Table_blastp_Pf3D7.iterrows():
    Table_blastp_Pf3D7.at[index, 'Absolute start'] = int(row['q. start']) - (int(row['s. start']) - 1)
    Table_blastp_Pf3D7.at[index, 'Absolute end'] = int(row['q. end']) + (11 - int(row['s. end']))
    if row['subject id'] in mhc_dict.keys():
        Table_blastp_Pf3D7.at[index, 'HLA recognition'] = mhc_dict[row['subject id']]
    else:
        Table_blastp_Pf3D7.at[index, 'HLA recognition'] = 'Not recognise by TopABC-3'
    if (Table_blastp_Pf3D7.at[index, 'Absolute end'] in netchop_dict.keys() and
            (netchop_dict[Table_blastp_Pf3D7.at[index, 'Absolute end']] == dict_seq_C_terminal[row['subject id']])):
        Table_blastp_Pf3D7.at[index, 'C terminal chop'] = netchop_dict[Table_blastp_Pf3D7.at[index, 'Absolute end']]
    else:
        Table_blastp_Pf3D7.at[index, 'C terminal chop'] = 'No'

FinalTable = Table_blastp_Pf3D7.copy()
FinalTable = FinalTable.drop(columns=['query id', 'bit score', 'evalue'])

# FinalTable.to_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
#                  'Cterminalmatches_location.csv')
