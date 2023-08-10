import pandas as pd
import pickle

from tqdm import tqdm

# Import of data and filtering based on rank

mhc_run = pd.read_table('/Users/u2176312/OneDrive - University of '
                        'Warwick/CSP/AllelePops/Filtered_HLAs_all_all_lenght_corrected.txt', sep='\t')
temp_file = open('/Users/u2176312/OneDrive - University of '
                 'Warwick/CSP/AllelePops/Kmer_CSP_region_273-375_aa.fasta_dict.pickle', 'rb')
mhc_run_dict = pickle.load(temp_file)

mhc_run_successful = mhc_run[mhc_run['rank'] <= 1]
mhc_run_successful = mhc_run_successful.sort_values(by=['allele'])

successful_alleles = list(set(mhc_run_successful['allele']))
successful_alleles.sort()

# %%
dict_ids = {}
starting_id = 0
final_dict = {}

for key in tqdm(mhc_run_dict.keys()):
    value = mhc_run_dict[key]
    if isinstance(value, str):
        number_of_sequences = 1
    else:
        number_of_sequences = len(mhc_run_dict[key])
    dict_ids.update({key: list(range(starting_id, starting_id + number_of_sequences))})
    starting_id = starting_id + number_of_sequences
    final_dict.update({key: number_of_sequences})

temp_dict = {}
for allele in tqdm(successful_alleles):
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['seq_num'] not in temp_dict.keys():
            temp_dict.update({row['seq_num']: [row['allele']]})

        else:
            temp_dict[row['seq_num']].append(row['allele'])

final_df = pd.DataFrame.from_dict(final_dict, orient='index', columns=['Variants'])

#%%
newdata = {col: [] * len(final_df) for col in successful_alleles}
new_df = pd.DataFrame(newdata)

updateddf = pd.concat([final_df, new_df], axis=1)

final_df = updateddf.copy()
final_df = final_df.fillna(0)

for key1 in tqdm(dict_ids.keys()):
    values = dict_ids[key1]
    for value in values:
        if value in temp_dict.keys():
            alleles = temp_dict[value]
            for i in alleles:
                final_df.loc[key1, i] = final_df.loc[key1, i] + 1
        else:
            continue

# %%
dict_tosave = {}
for allele in successful_alleles:
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['allele'] not in dict_tosave.keys():
            dict_tosave.update({row['allele']: [int(row['seq_num']) - 1]})

        else:
            dict_tosave[row['allele']].append(int(row['seq_num']) - 1)

sheet2 = pd.DataFrame.from_dict(dict_tosave, orient='index')
sheet2 = sheet2.transpose()

writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of Warwick/CSP/Kmer_CSP_region_283-307_summarydata.xlsx',
                        engine='openpyxl')
final_df.to_excel(writer, sheet_name='summarydata')
sheet2.to_excel(writer, sheet_name='Alleles&Sequences')
writer.close()

with open('/Users/u2176312/OneDrive - University of Warwick/CSP/Allele_seq_ids.pickle', 'wb') as f:
    pickle.dump(dict_tosave, f)


# %% Dictionary of the NCBI sequences and kmers
ncbi_prediction = pd.read_table('/Users/u2176312/OneDrive - University of '
                                'Warwick/CSP/NCBI_CSP/TopABC_matchingseqs/TopABC_matchingseqs.txt')
mhc_run_successful = ncbi_prediction.sort_values(by=['allele'])

successful_alleles = list(set(mhc_run_successful['allele']))
successful_alleles.sort()
temp_file = open('/Users/u2176312/OneDrive - University of '
                 'Warwick/CSP/NCBI_CSP/NCBI_CSP_peptides_11kmer.pickle', 'rb')
ncbi_dict = pickle.load(temp_file)

dict_ids = {}
starting_id = 0
final_dict = {}
kmer_start = 0
kmer_end = 11

for key in ncbi_dict.keys():
    value = ncbi_dict[key]
    if isinstance(value, str):
        number_of_sequences = 1
    else:
        number_of_sequences = len(ncbi_dict[key])

    kmerid = 'Kmer_' + str(kmer_start) + '_' + str(kmer_end)
    dict_ids.update({kmerid: list(range(starting_id, starting_id + number_of_sequences))})
    starting_id = starting_id + number_of_sequences
    final_dict.update({kmerid: number_of_sequences})
    kmer_start += 1
    kmer_end += 1


temp_dict = {}
for allele in successful_alleles:
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['seq_num'] not in temp_dict.keys():
            temp_dict.update({row['seq_num']: [row['allele']]})

        else:
            temp_dict[row['seq_num']].append(row['allele'])

final_df = pd.DataFrame.from_dict(final_dict, orient='index', columns=['Variants'])

for allele in successful_alleles:
    final_df[allele] = 0

for key1 in dict_ids.keys():
    values = dict_ids[key1]
    for value in values:
        if value in temp_dict.keys():
            alleles = temp_dict[value]
            for i in alleles:
                final_df.loc[key1, i] = final_df.loc[key1, i] + 1
        else:
            continue

dict_tosave = {}
for allele in successful_alleles:
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['allele'] not in dict_tosave.keys():
            dict_tosave.update({row['allele']: [int(row['seq_num']) - 1]})

        else:
            dict_tosave[row['allele']].append(int(row['seq_num']) - 1)

sheet2 = pd.DataFrame.from_dict(dict_tosave, orient='index')
sheet2 = sheet2.transpose()

writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
                        'TopABC_matchingseqs/TopABC_matchingseqs_summarydata.xlsx',
                        engine='openpyxl')
final_df.to_excel(writer, sheet_name='summarydata')
sheet2.to_excel(writer, sheet_name='Alleles&Sequences')
writer.close()

with open('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
          'TopABC_matchingseqs/TopABC_matchingseqs.pickle', 'wb') as f:
    pickle.dump(dict_tosave, f)
