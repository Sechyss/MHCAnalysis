import pandas as pd
import pickle

mhc_kenya_run = pd.read_table('/Users/u2176312/OneDrive - University of '
                              'Warwick/CSP/CSP_SNP_region/Kmer_CSP_predictionbinding_region_308-342.txt')
temp_file = open('/Users/u2176312/OneDrive - University of '
                 'Warwick/CSP/CSP_SNP_region/Kmer_CSP_region_283-307_aa.fasta_dict.pickle', 'rb')
mhc_kenya_run_dict = pickle.load(temp_file)

mhc_kenya_run_successful = mhc_kenya_run[mhc_kenya_run['rank'] <= 1]
mhc_kenya_run_successful = mhc_kenya_run_successful.sort_values(by=['allele'])

successful_alleles = list(set(mhc_kenya_run_successful['allele']))
successful_alleles.sort()

dict_ids = {}
starting_id = 0
final_dict = {}

for key in mhc_kenya_run_dict.keys():
    number_of_sequences = len(mhc_kenya_run_dict[key])
    dict_ids.update({key: list(range(starting_id, starting_id + number_of_sequences))})
    starting_id = starting_id + number_of_sequences
    final_dict.update({key: number_of_sequences})

temp_dict = {}
for allele in successful_alleles:
    filtered_df = mhc_kenya_run_successful[mhc_kenya_run_successful['allele'] == allele]
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
                final_df.loc[key1][i] += 1
        else:
            continue

dict_tosave = {}
for allele in successful_alleles:
    filtered_df = mhc_kenya_run_successful[mhc_kenya_run_successful['allele'] == allele]
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
