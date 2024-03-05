import pandas as pd
import pickle

from tqdm import tqdm
from Bio import SeqIO


def slice_dict(dictionary, keys):
    """
  Slices a dictionary to extract only keys and values that match a coordinates of keys.

  Args:
    dictionary: The dictionary to slice.
    keys: The coordinates of keys to match.

  Returns:
    A new dictionary containing only the keys and values that match the coordinates of keys.
  """
    new_dict = {}
    for key2 in keys:
        if key2 in dictionary:
            new_dict[key2] = dictionary[key2]
    return new_dict


# %%  Import of data and filtering based on rank

mhc_run = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                        'CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/NCBI_TopABC_lenght11_ncbiseqs.txt',
                        sep='\t')
temp_file = open('/Users/u2176312/OneDrive - University of '
                 'Warwick/CSP/NCBI_CSP/NCBI_CSP_peptides_11kmer.pickle', 'rb')
mhc_run_dict = pickle.load(temp_file)

mhc_run_successful = mhc_run[mhc_run['rank'] <= 1]
mhc_run_successful = mhc_run_successful.sort_values(by=['allele'])

successful_alleles = list(set(mhc_run_successful['allele']))
successful_alleles.sort()

# %% Dicitonary of sequence order

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                        'CSP/NCBI_CSP/NCBI_CSP_peptides_11kmer_filtered.fasta', 'fasta')

dict_sequence_order = {}
for seq_record in fastafile:
    sequence = str(seq_record.seq)
    number = int(seq_record.id.replace('Sequence_', ''))
    dict_sequence_order.update({sequence: number})

# %%

final_dict = {}

for key1 in tqdm(mhc_run_dict.keys()):
    value = mhc_run_dict[key1]
    if isinstance(value, str):
        number_of_sequences = 1
    else:
        number_of_sequences = len(mhc_run_dict[key1])
    final_dict.update({key1: number_of_sequences})

matching_dict = slice_dict(final_dict, dict_sequence_order.keys())

temp_dict = {}
for allele in tqdm(successful_alleles):
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['seq_num'] not in temp_dict.keys():
            temp_dict.update({row['seq_num']: [row['allele']]})

        else:
            temp_dict[row['seq_num']].append(row['allele'])

final_df = pd.DataFrame.from_dict(matching_dict, orient='index', columns=['Number of genomes'])

# %%
newdata = {col: [] * len(final_df) for col in successful_alleles}
new_df = pd.DataFrame(newdata)

updateddf = pd.concat([final_df, new_df], axis=1)

final_df = updateddf.copy()
final_df['HLAs recognised'] = ''
final_df = final_df.fillna(0)

for key in dict_sequence_order.keys():
    sequence_number = dict_sequence_order[key]
    if sequence_number in temp_dict.keys():
        list_hlas = list(set(temp_dict[sequence_number]))

        for hla in list_hlas:
            final_df.loc[key, hla] = final_df.loc[key, hla] + 1
final_df['HLAs recognised'] = final_df.loc[:, final_df.columns.str.contains('HLA-')].sum(axis=1)

# %%
dict_tosave = {}
for allele in successful_alleles:
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['allele'] not in dict_tosave.keys():
            dict_tosave.update({row['allele']: [str(int(row['seq_num']) - 1) + '_length_' + str(row['length'])]})

        else:
            dict_tosave[row['allele']].append(str(int(row['seq_num']) - 1) + '_length_' + str(row['length']))

sheet2 = pd.DataFrame.from_dict(dict_tosave, orient='index')
sheet2 = sheet2.transpose()

writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
                        'resultsPredictionBinding_NCBI_only_TopABC/TopABC_NCBI_length11_summarydata.xlsx',
                        engine='openpyxl')
final_df.to_excel(writer, sheet_name='summarydata')
sheet2.to_excel(writer, sheet_name='Alleles&Sequences')
writer.close()

with open('/Users/u2176312/OneDrive - University of Warwick/'
          'CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/Allele_seq_ids_NCBI_length11.pickle', 'wb') as f:
    pickle.dump(temp_dict, f)
