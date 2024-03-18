# /usr/bin/env python3
import pickle

import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

# Define column names for the blastp result table
columns = ['query id', 'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens',
           'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']

# %%  Upload of the different datasets that will be used later in the code
#  Blastp result table related to human recognition only focusing those not similar to human peptides
Blastp_human_recognition = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                                         'Otherproteins/BLASTP_results/Blastp_Human_CSP_Nterminal.tsv', header=None)
Blastp_human_recognition.columns = columns

# Modify the 'subject id' column by removing 'Sequence_' and adding 1 to the values
Blastp_human_recognition['subject id'] = Blastp_human_recognition['subject id'].apply(
    lambda x: int(str(x).replace('Sequence_', '')) + 1)
dictionary_human = Blastp_human_recognition.groupby("subject id")['% identity'].apply(set).to_dict()
Blastp_human_recognition = Blastp_human_recognition[Blastp_human_recognition['% identity'] >= 80]
HumanKmers = list(set(Blastp_human_recognition['subject id']))

# Read a table of MHC data from a TSV file
Table_mhc = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                          'Otherproteins/IEDB_prediction/NCBI_CSP_Nterminal_IEDB.tsv', sep='\t')
Table_mhc.dropna(how='all', inplace=True)
Table_mhc['seq_num'] = Table_mhc['seq_num'].astype(int)

mhc_dict = {}
for index, row in Table_mhc.iterrows():
    if row['seq_num'] not in mhc_dict.keys():
        mhc_dict.update({row['seq_num']: [[row['allele']], [row['start']], [row['end']],
                                          [row['peptide']], [str(row['peptide'])[-1]]]})
    else:
        mhc_dict[row['seq_num']][0].append(row['allele'])
        mhc_dict[row['seq_num']][1].append(row['start'])
        mhc_dict[row['seq_num']][2].append(row['end'])
        mhc_dict[row['seq_num']][3].append(row['peptide'])
        mhc_dict[row['seq_num']][4].append(str(row['peptide'])[-1])

# Creation of sequence dictionary

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                        'Otherproteins/Multiple_Alignment/CSPNterminal_kmer_filtered_corrected.fasta',
                        'fasta')

dictionary_sequences = {}

for seq_record in fastafile:
    sequence_id = int(str(seq_record.id).replace('Sequence_', '')) + 1
    peptide_sequence = seq_record.seq
    dictionary_sequences.update({sequence_id: peptide_sequence})

dictionary_depth = pickle.load(open('/Users/u2176312/OneDrive - University of Warwick/'
                               'Otherproteins/Multiple_Alignment/CSPNterminal_Kmer_depth.pickle', 'rb'))

dictionary_locations = pickle.load(open('/Users/u2176312/OneDrive - University of Warwick/'
                                        'Otherproteins/Multiple_Alignment/CSPNterminal.pickle', 'rb'))
Table_blastp_Pf3D7 = pd.DataFrame(columns=['Seq_id', 'Absolute start', 'Absolute end'
    , 'Peptide of 11kmer-aa'])
counter_1 = 1
for key in dictionary_sequences.keys():
    sequence_kmer = dictionary_sequences[key]
    counter_2 = 0
    for values in dictionary_locations[sequence_kmer]:
        if int(str(values).rsplit('_')[1]) != counter_2:
            location_start = int(str(values).rsplit('_')[1])
            location_end = int(str(values).rsplit('_')[2])
            Table_blastp_Pf3D7.loc[counter_1, 'Seq_id'] = key
            Table_blastp_Pf3D7.loc[counter_1, 'Absolute start'] = location_start
            Table_blastp_Pf3D7.loc[counter_1, 'Absolute end'] = location_end
            Table_blastp_Pf3D7.loc[counter_1, 'Peptide of 11kmer-aa'] = sequence_kmer

            counter_1 += 1
            counter_2 = location_start
df_without_duplicates = Table_blastp_Pf3D7.drop_duplicates()
# %% Creation of the first database using the previous tables

collectordf = pd.DataFrame(columns=['Sequence', 'Allele', 'Starting aa', 'Ending aa', 'Epitope', 'C-terminal aa'])

print('Analyzing the sequences in the dictionary mhc data')
for key in tqdm(mhc_dict.keys()):
    for element in range(len(mhc_dict[key][0])):
        row_to_add = pd.DataFrame({'Sequence': key,
                                   'Allele': mhc_dict[key][0][element],
                                   'Starting aa': mhc_dict[key][1][element],
                                   'Ending aa': mhc_dict[key][2][element],
                                   'Epitope': mhc_dict[key][3][element],
                                   'C-terminal aa': mhc_dict[key][4][element]}, index=[0])
        collectordf = pd.concat([row_to_add, collectordf.loc[:]]).reset_index(drop=True)

# Preparation of the final dataframe which contains all the relevant data

final_df = collectordf.set_index('Sequence').join(df_without_duplicates.set_index('Seq_id'), how='outer')
final_df.reset_index(inplace=True)
final_df = final_df.dropna(how='all')

final_df['Human peptide recognition'] = ''
final_df['% identity with human peptidome'] = ''
for index, row in final_df.iterrows():

    if row['Sequence'] in dictionary_human.keys():
        final_df.at[index, '% identity with human peptidome'] = float(max(dictionary_human[row['Sequence']]))
    else:
        final_df.at[index, '% identity with human peptidome'] = 0

    if row['Sequence'] in HumanKmers:
        final_df.at[index, 'Human peptide recognition'] = 1
    else:
        final_df.at[index, 'Human peptide recognition'] = 0


final_df.dropna(subset=['Absolute start'], inplace=True)
Filtered_df = final_df[final_df['Human peptide recognition'] == 0]
Table = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/'
                      'CSP/AllelePops/FilteredDataAllele.xlsx', sheet_name='TOP_ABC')

result_dict = Table.to_dict(orient='dict')

list_countries = [x for x in result_dict.keys()]
color_codes = [
    '#E69F00',
    '#56B4E9',
    '#009E73',
    '#F0E442',
    '#0072B2',
    '#D55E00',
    '#CC79A7',
    '#88CCEE',
    '#CC6677',
    '#DDCC77',
    '#117733',
    '#332288',
    '#AA4499',
    '#44AA99',
    '#999933',
    '#882255',
    '#661100',
    '#6699CC',
    '#888888',
]
color_dictionary = dict(zip(list_countries, color_codes))

new_df2 = pd.DataFrame()
print('Analysing the different countries and collecting the data')
for country in tqdm(result_dict.keys()):
    x_axis = []
    number_variants = []
    hla_recognition = []
    sequences_recognition = []
    number_variations = []
    min_similarity_hum_peptidome = []
    max_similarity_hum_peptidome = []
    average_similarity_hum_peptidome = []

    listHLAs = list(result_dict[country].values())
    countryDF = Filtered_df[Filtered_df['Allele'].isin(listHLAs)].sort_values(by=['Absolute start'], ascending=True)
    startingkmer = final_df['Absolute start'].min()
    endingkmer = startingkmer + 10

    while startingkmer <= final_df['Absolute start'].max():
        x_axis.append('Kmer_' + str(int(startingkmer)) + '_' + str(int(endingkmer)))  # Creation the kmer axis

        tempdf = countryDF[countryDF['Absolute start'] == startingkmer]  # Filtering the country df to the kmer axis
        temp_list = tempdf['Absolute start'].tolist()
        if len(temp_list) > 0:
            alleles = list(set(tempdf['Allele'].to_list()))  # List of HLAs in the countrydf
            sequences_country = set(tempdf['Sequence'].to_list())  # List of variation in the selected countrydf
            hla_recognition.append(len(alleles))
            sequences_recognition.append(len(sequences_country))
        else:
            hla_recognition.append(0)
            sequences_recognition.append(0)

        tempdf2 = Filtered_df[Filtered_df['Absolute start'] == startingkmer]  # selection of the filtered df to kmer
        temp_list = tempdf2['Absolute start'].tolist()
        if len(temp_list) > 0:
            sequences = set(tempdf2['Sequence'].to_list())
            number_variants.append(len(sequences))
        else:
            number_variants.append(0)

        tempdf3 = final_df[final_df['Absolute start'] == startingkmer]
        temp_list = tempdf3['Absolute start'].tolist()
        if len(temp_list) > 0:
            variations = set(tempdf3['Sequence'].astype(str).to_list())
            min_similarity_hum_peptidome.append(min(set(tempdf3['% identity with human peptidome'].to_list())))
            max_similarity_hum_peptidome.append(max(set(tempdf3['% identity with human peptidome'].to_list())))
            average_similarity_hum_peptidome.append(
                (sum(set(tempdf3['% identity with human peptidome'].to_list())) /
                 len(set(tempdf3['% identity with human peptidome'].to_list()))))

            number_variations.append(len(variations))

        else:
            min_similarity_hum_peptidome.append(0)
            max_similarity_hum_peptidome.append(0)
            average_similarity_hum_peptidome.append(0)

            number_variations.append(0)

        startingkmer += 1
        endingkmer += 1

    new_df2.index = x_axis

    new_df2['Total variations'] = number_variations
    new_df2['Coverage depth'] = list(dictionary_depth.values())
    new_df2['Kmer min% human peptidome'] = min_similarity_hum_peptidome
    new_df2['Kmer average% human peptidome'] = average_similarity_hum_peptidome
    new_df2['Kmer max% human peptidome'] = max_similarity_hum_peptidome
    new_df2[country] = sequences_recognition
