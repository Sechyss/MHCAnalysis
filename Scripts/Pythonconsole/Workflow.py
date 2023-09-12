# Import necessary libraries
import pandas as pd  # Pandas library for data manipulation
import pickle
from dna_features_viewer import GraphicFeature, GraphicRecord
from tqdm import tqdm

# Define column names for the blastp result table
columns = ['query id', 'subject id', '% identity', 'alignment length', 'mismatches', 'gap opens',
           'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score']

#  Blastp result table related to human recognition only focusing those not similar to human peptides
Blastp_human_recognition = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/'
                                         'CSP/Humanrecogn/Human_Malaria_NCBI_blastp.tsv', header=None)
Blastp_human_recognition.columns = columns
Blastp_human_recognition = Blastp_human_recognition[Blastp_human_recognition['% identity'] >= 80]
# Modify the 'subject id' column by removing 'Sequence_' and adding 1 to the values
Blastp_human_recognition['subject id'] = Blastp_human_recognition['subject id'].apply(
    lambda x: int(str(x).replace('Sequence_', '')) + 1)
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

# Remove duplicate values within the MHC data lists
for key in mhc_dict.keys():
    new_data = list(set(mhc_dict[key][0]))
    mhc_dict[key][0] = new_data

# Read blastp result table from a TSV file into a Pandas DataFrame
Table_blastp_Pf3D7 = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                                   'NCBI_CSP/NCBI_Pf3D7_blastp.tsv', sep='\t', header=None)
Table_blastp_Pf3D7.columns = columns  # Assign column names

# Modify the 'subject id' column by removing 'Sequence_' and adding 1 to the values
Table_blastp_Pf3D7['subject id'] = Table_blastp_Pf3D7['subject id'].apply(
    lambda x: int(str(x).replace('Sequence_', '')) + 1)

Table_blastp_Pf3D7 = Table_blastp_Pf3D7.drop(['query id'], axis=1)

# Read Netchop result table from a CSV file
NetchopTable = pd.read_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/Netchop/result_CSP_Pf3D7_netchop.csv')
NetchopTable = NetchopTable[NetchopTable['prediction_score'] >= 0.5]  # Filter rows with prediction score >= 0.5

# Create a dictionary to store Netchop data for each '#'
netchop_dict = {}
for index, row in NetchopTable.iterrows():
    netchop_dict.update({row['#']: row['amino_acid']})

# %% Collection of data from the dictionary for MHC data
collectordf = pd.DataFrame(columns=['Sequence', 'Allele', 'Starting aa', 'Ending aa', 'Peptide', 'C-terminal aa'])

for key in tqdm(mhc_dict.keys()):
    for element in range(len(mhc_dict[key][0])):
        row_to_add = pd.DataFrame({'Sequence': key,
                                   'Allele': mhc_dict[key][0][element],
                                   'Starting aa': mhc_dict[key][1][element],
                                   'Ending aa': mhc_dict[key][2][element],
                                   'Peptide': mhc_dict[key][3][element],
                                   'C-terminal aa': mhc_dict[key][4][element]}, index=[0])
        collectordf = pd.concat([row_to_add, collectordf.loc[:]]).reset_index(drop=True)

# %% Preparation of the final dataframe which contains all the relevant data

final_df = collectordf.set_index('Sequence').join(Table_blastp_Pf3D7.set_index('subject id'), how='inner')
final_df = final_df.drop(['alignment length', 'bit score', 'gap opens', 'evalue'], axis=1)

final_df['Human peptide recognition'] = ''
final_df['Absolute end (MHC peptide)'] = ''
final_df['C-terminal match'] = ''
for index, row in final_df.iterrows():
    end_MHC_peptide = row['Ending aa']
    end_blastpPf3D7 = row['s. end']
    absolute_end = row['q. end'] - (end_blastpPf3D7 - end_MHC_peptide)
    final_df.at[index, 'Absolute end (MHC peptide)'] = absolute_end
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

final_df.to_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
                'Cterminalmatches_location.csv', index_label='Sequence number')

# %% Plotting of the results
final_df = final_df[(final_df['C-terminal match'] == 1) & (final_df['Human peptide recognition'] == 1)]

tempfile = open('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Dictionary_country_alleles.pickle',
                'rb')
dictionaryHLA_countries = pickle.load(tempfile)

# Features includes the name of the sequence, the length of the sequence, the number of elements...
features = [GraphicFeature(start=1, end=18, color="olivedrab", label='Signal'),
            GraphicFeature(start=93, end=97, color="brown", label='Region I'),
            GraphicFeature(start=105, end=272, color="skyblue", label='Tandem repeats'),
            GraphicFeature(start=326, end=342, color="crimson", label='Th2R region'),
            GraphicFeature(start=366, end=380, color='green', label='Th3R region'),
            GraphicFeature(start=375, end=397, color="purple", label='GPI-anchor')]

for country in dictionaryHLA_countries.keys():
    startingkmer = 0
    endingKmer = 0
    listHLAs = dictionaryHLA_countries[country]
    listHLAs = ['HLA-'+str(x) for x in listHLAs]
    countryDF = final_df[final_df['Allele'].isin(listHLAs)].sort_values(by=['q. start'], ascending=True)
    for index, row in countryDF.iterrows():
        if row['q. start'] > endingKmer:
            startingkmer = row['q. start']
            endingKmer = startingkmer + 11
            tempdf = countryDF[(countryDF['q. start'] >= startingkmer) &
                               (countryDF['Absolute end (MHC peptide)'] <= endingKmer)]
            alleles = list(set(tempdf['Allele'].to_list))
            features.append(GraphicFeature(start=startingkmer, end=endingKmer, color='grey', label=str(len(alleles))))



