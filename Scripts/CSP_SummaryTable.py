#!/usr/bin/env python3

# Upload of packages to run the script
import argparse
import pickle
import textwrap
import warnings


import pandas as pd


warnings.simplefilter(action='ignore', category=FutureWarning)
# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Summary table of the prediction binding protein in combination with the Sequences generated by the combination of all 
the SNPs in fasta file.
------------------------------------------
'''))

parser.add_argument("-t", "--table", metavar='table.txt', dest="mhctable", help="MHC prediction binding table",
                    type=str)
parser.add_argument("-d", "--dictionary", metavar='dict.pickle', dest="dictionary", help="Pickle dictionary", type=str)
parser.add_argument("-p", "--newpickle", dest="newdict", help="Allele Id pickle dictionary", type=str)
parser.add_argument("-o", "--output", metavar='file.xlsx', dest="output", help="Output file to save", type=str)

args = parser.parse_args()

# =============================================================================

# Upload the information from table and the dictionary
mhc_run = pd.read_table(args.mhctable)
temp_file = open(args.dictionary, 'rb')
mhc_run_dict = pickle.load(temp_file)

# Filter results to only those with a binding site below 1 and sort them by allele id
mhc_run_successful = mhc_run[mhc_run['rank'] <= 1]
mhc_run_successful = mhc_run_successful.sort_values(by=['allele'])

# Extract the list of alleles
successful_alleles = list(set(mhc_run_successful['allele']))
successful_alleles.sort()

# Creation of constants to be used in the script
dict_ids = {}
starting_id = 0
final_dict = {}

# Association of Kmer with the sequences ids
for key in mhc_run_dict.keys():
    number_of_sequences = len(mhc_run_dict[key])
    dict_ids.update({key: list(range(starting_id, starting_id + number_of_sequences))})
    starting_id = starting_id + number_of_sequences
    final_dict.update({key: number_of_sequences})

temp_dict = {}
# Link between Sequences and alleles that bound successfully
for allele in successful_alleles:
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['seq_num'] not in temp_dict.keys():
            temp_dict.update({row['seq_num']: [row['allele']]})

        else:
            temp_dict[row['seq_num']].append(row['allele'])

# Creation of the final dataframe
final_df = pd.DataFrame.from_dict(final_dict, orient='index', columns=['Variants'])

for allele in successful_alleles:
    final_df[allele] = 0

# Sum the number of sequences in each allele that bind to a particular Kmer
for key1 in dict_ids.keys():
    values = dict_ids[key1]
    for value in values:
        if value in temp_dict.keys():
            alleles = temp_dict[value]
            for i in alleles:
                final_df.loc[key1][i] += 1
        else:
            continue

# Creation of dictionary to save the list of alleles and the sequence ids they bound to
dict_tosave = {}
for allele in successful_alleles:
    filtered_df = mhc_run_successful[mhc_run_successful['allele'] == allele]
    redefined_df = filtered_df.sort_values(by=['seq_num'])
    for index, row in redefined_df.iterrows():
        if row['allele'] not in dict_tosave.keys():
            dict_tosave.update({row['allele']: [int(row['seq_num'])-1]})

        else:
            dict_tosave[row['allele']].append(int(row['seq_num'])-1)

sheet2 = pd.DataFrame.from_dict(dict_tosave, orient='index')
sheet2 = sheet2.transpose()

# Creation of the Excel spreadsheet to save the table
writer = pd.ExcelWriter(args.output, engine='openpyxl')
final_df.to_excel(writer, sheet_name='summarydata')
sheet2.to_excel(writer, sheet_name='Alleles&Sequences')
writer.close()

# Save the dictionary as pickle file
with open(args.newdict, 'wb') as f:
    pickle.dump(dict_tosave, f)
