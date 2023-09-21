# Import necessary libraries
import pandas as pd  # Pandas library for data manipulation

from matplotlib import pyplot as plt
from tqdm import tqdm
from Bio import SeqIO

# Define column names for the blastp result table
columns = ['query id', 'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens',
           'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']

#  Blastp result table related to human recognition only focusing those not similar to human peptides
Blastp_human_recognition = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                                         'CSP/Humanrecogn/Human_Malaria_NCBI_blastp.tsv', header=None)
Blastp_human_recognition.columns = columns

# Modify the 'subject id' column by removing 'Sequence_' and adding 1 to the values
Blastp_human_recognition['subject id'] = Blastp_human_recognition['subject id'].apply(
    lambda x: int(str(x).replace('Sequence_', '')) + 1)
dictionary_human = Blastp_human_recognition.groupby("subject id")['% identity'].apply(set).to_dict()
Blastp_human_recognition = Blastp_human_recognition[Blastp_human_recognition['% identity'] >= 80]
HumanKmers = list(set(Blastp_human_recognition['subject id']))

# Read a table of MHC data from a TSV file and filter it to include only rows with 'seq_num' in Highpercentage
Table_mhc = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                          'NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
                          'NCBI_TopABC_all_lengths_NCBIseqs.txt', sep='\t')

# Create a dictionary to store MHC data for each 'seq_num'
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

# Read blastp result table from a TSV file into a Pandas DataFrame
Table_blastp_Pf3D7 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                                   'NCBI_CSP/Blastp_NCBI_Pf3D7.tsv', sep='\t', header=None)
Table_blastp_Pf3D7.columns = columns  # Assign column names
Table_blastp_Pf3D7 = Table_blastp_Pf3D7[Table_blastp_Pf3D7['s. start'] > 270]

# Modify the 'subject id' column by removing 'Sequence_' and adding 1 to the values
Table_blastp_Pf3D7['query id'] = Table_blastp_Pf3D7['query id'].apply(
    lambda x: int(str(x).replace('Sequence_', '')) + 1)

Table_blastp_Pf3D7 = Table_blastp_Pf3D7.drop(['subject id', 'alignment length',
                                              'bit score', 'gap opens', 'evalue'], axis=1)

# Read Netchop result table from a CSV file
NetchopTable = pd.read_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/Netchop/result_CSP_Pf3D7_netchop.csv')
NetchopTable = NetchopTable[NetchopTable['prediction_score'] >= 0.5]  # Filter rows with prediction score >= 0.5

# Create a dictionary to store Netchop data for each '#'
netchop_dict = {}
for index, row in NetchopTable.iterrows():
    netchop_dict.update({row['#']: row['amino_acid']})

# %% Creation of sequence dictionary

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
                        'NCBI_CSP_peptides_11kmer_filtered.fasta', 'fasta')

dictionary_sequences = {}

for seq_record in fastafile:
    sequence_id = int(str(seq_record.id).replace('Sequence_', '')) + 1
    peptide_sequence = seq_record.seq
    dictionary_sequences.update({sequence_id: peptide_sequence})

# %% Collection of data from the dictionary for MHC data
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

# %% Preparation of the final dataframe which contains all the relevant data

final_df = collectordf.set_index('Sequence').join(Table_blastp_Pf3D7.set_index('query id'), how='inner')
final_df.reset_index(inplace=True)

final_df['Peptide of 11kmer-aa'] = ''
final_df['Absolute end (MHC peptide)'] = ''
final_df['C-terminal match'] = ''
final_df['Human peptide recognition'] = ''
final_df['% identity with human peptidome'] = ''
for index, row in final_df.iterrows():
    end_MHC_peptide = row['Ending aa']
    end_blastpPf3D7 = row['q. end']
    absolute_end = row['s. end'] - (end_blastpPf3D7 - end_MHC_peptide)
    final_df.at[index, 'Absolute end (MHC peptide)'] = absolute_end
    if row['Sequence'] in dictionary_human.keys():
        final_df.at[index, '% identity with human peptidome'] = str(max(dictionary_human[row['Sequence']]))
    else:
        final_df.at[index, '% identity with human peptidome'] = 0
    final_df.at[index, 'Peptide of 11kmer-aa'] = str(dictionary_sequences[row['Sequence']])
    if absolute_end in netchop_dict.keys():
        C_value = netchop_dict[absolute_end]
        if C_value == row['C-terminal aa']:
            final_df.at[index, 'C-terminal match'] = 1
        else:
            final_df.at[index, 'C-terminal match'] = 0
    else:
        final_df.at[index, 'C-terminal match'] = 'No Netchop at Absolute end position'

    if index in HumanKmers:
        final_df.at[index, 'Human peptide recognition'] = 1
    else:
        final_df.at[index, 'Human peptide recognition'] = 0

final_df.to_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/AllelePopNCBI_Workflow/'
                'Cterminalmatches_location.csv', index=False)

# %% Plotting of the results
Filtered_df = final_df[(final_df['C-terminal match'] == 1) & (final_df['Human peptide recognition'] == 0)]

Table_blastp_Pf3D7['Absolute start'] = ''
Table_blastp_Pf3D7['Absolute end'] = ''

for index, row in Table_blastp_Pf3D7.iterrows():
    Table_blastp_Pf3D7.at[index, 'Absolute start'] = row['s. start'] - (row['q. start'] - 1)
    Table_blastp_Pf3D7.at[index, 'Absolute end'] = row['s. end'] - (row['q. end'] - 11)

Table = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/FilteredDataAllele.xlsx',
                      sheet_name='TOP_ABC')

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

new_df = pd.DataFrame()
new_df2 = pd.DataFrame()
print('Generating the figures for the different countries')
for country in tqdm(result_dict.keys()):
    x_axis = []
    number_variants = []
    hla_recognition = []
    sequences_recognition = []
    number_variations = []

    listHLAs = list(result_dict[country].values())
    countryDF = Filtered_df[Filtered_df['Allele'].isin(listHLAs)].sort_values(by=['s. start'], ascending=True)
    startingkmer = Filtered_df['s. start'].min()
    endingkmer = startingkmer + 11
    while endingkmer <= Filtered_df['Absolute end (MHC peptide)'].max():
        x_axis.append('Kmer_' + str(startingkmer) + '_' + str(endingkmer))
        tempdf = countryDF[(countryDF['s. start'] >= startingkmer) &
                           (countryDF['Absolute end (MHC peptide)'] <= endingkmer)]
        alleles = list(set(tempdf['Allele'].to_list()))
        sequences_country = list(set(tempdf['Sequence'].to_list()))
        hla_recognition.append(len(alleles))
        tempdf2 = Filtered_df[(Filtered_df['s. start'] >= startingkmer) &
                              (Filtered_df['Absolute end (MHC peptide)'] <= endingkmer)].sort_values(by=['s. start'],
                                                                                                     ascending=True)
        sequences = list(set(tempdf2['Sequence'].to_list()))
        number_variants.append(len(sequences))

        tempdf3 = Table_blastp_Pf3D7[(Table_blastp_Pf3D7['Absolute start'] >= startingkmer) &
                                     (Table_blastp_Pf3D7['Absolute end'] <= endingkmer)]
        variations = list(set(tempdf3['query id'].to_list()))
        variations = list(set(variations).union(set(sequences)))
        sequences_recognition.append(len(sequences_country))
        number_variations.append(len(variations))

        startingkmer += 1
        endingkmer += 1
    new_df.index = x_axis
    new_df2.index = x_axis
    new_df['Sequences with C-terminal / NonHuman'] = number_variants
    new_df2['Total variations'] = number_variations
    new_df[country] = hla_recognition
    new_df2[country] = sequences_recognition
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(25, 20), sharex=True)
    # Creation of the plot with the different populations
    ax1.bar(x=x_axis, height=number_variations, edgecolor='black', color='white')
    ax1.set_title('Number of variations per Kmer', loc='left', y=1, weight='bold', fontsize=30)
    ax2.bar(x=x_axis, height=sequences_recognition, edgecolor='black', color=color_dictionary[country])
    ax2.set_title(country + ' number of variants recognised by HLAs', loc='left', y=1, weight='bold', fontsize=30)

    plt.xticks(fontsize=10, fontweight='bold', rotation='vertical')

    plt.tight_layout()
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/AllelePopNCBI_Workflow/'
                'Kmer_NCBI_workflow_recognition_' + country + '.pdf', dpi=600)


writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/AllelePopNCBI_Workflow/'
                        'AllelePopNCBI_WorkflowSummaryTable.xlsx', engine='openpyxl')
new_df.to_excel(writer, sheet_name='HLA numbers')
new_df2.to_excel(writer, sheet_name='RelativeFreq')

writer.close()
