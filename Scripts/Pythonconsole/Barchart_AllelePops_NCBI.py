# %% Import the different packages
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import pickle
import seaborn as sns


def add_labels(x, z, axis):  # Add labels to the plot later
    for element in range(len(x)):
        axis.text(element, z[element], z[element], ha='center', fontweight='bold', color='black', fontsize=5)



# %% List of TopABC but only the top 2, in case you want to filter them without repeating analysis
list_topABC_2 = [
    'HLA-A*01:01',
    'HLA-A*02:01',
    'HLA-A*02:06',
    'HLA-A*02:11',
    'HLA-A*03:01',
    'HLA-A*11:01',
    'HLA-A*23:01',
    'HLA-A*24:02',
    'HLA-A*30:01',
    'HLA-A*34:01',
    'HLA-B*07:02',
    'HLA-B*08:01',
    'HLA-B*13:01',
    'HLA-B*15:02',
    'HLA-B*15:06',
    'HLA-B*18:01',
    'HLA-B*35:01',
    'HLA-B*35:03',
    'HLA-B*40:01',
    'HLA-B*40:02',
    'HLA-B*40:06',
    'HLA-B*42:01',
    'HLA-B*44:02',
    'HLA-B*45:01',
    'HLA-B*46:01',
    'HLA-B*51:01',
    'HLA-B*52:01',
    'HLA-B*53:01',
    'HLA-B*54:01',
    'HLA-B*56:01',
    'HLA-B*56:02',
    'HLA-B*58:01',
    'HLA-C*01:02',
    'HLA-C*03:02',
    'HLA-C*03:03',
    'HLA-C*04:01',
    'HLA-C*04:03',
    'HLA-C*05:01',
    'HLA-C*06:02',
    'HLA-C*07:01',
    'HLA-C*07:02',
    'HLA-C*08:01',
    'HLA-C*08:01',
    'HLA-C*16:01',
    'HLA-C*16:01',
    'HLA-C*17:01',
]

# %% Upload the data and the functions to run the plot
Dataframe = pd.read_excel('/Users/u2176312/OneDrive - University of '
                          'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/TopABC_NCBI_summarydata.xlsx',
                          sheet_name='summarydata', index_col=0)
tempfile = open('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
                'resultsPredictionBinding_NCBI_only_TopABC/Allele_seq_ids_NCBI.pickle', 'rb')
dictionary_alleles = pickle.load(tempfile)
# dictionary_topABC_2 = {key: dictionary_alleles[key] for key in list_topABC_2}
x_axis = Dataframe.index.to_numpy()
number_genomes = Dataframe['Number of genomes'].to_numpy()
percentageoftopABC = Dataframe['HLAs recognised'].to_numpy()
percentageoftopABC = (percentageoftopABC / 56) * 100

# %%

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
# Creation of the plot with the different populations
ax1.bar(x=x_axis, height=number_genomes, edgecolor='black', color='white')
add_labels(x_axis, number_genomes, axis=ax1)
ax1.set_title('Number of genomes with peptides', loc='left', y=1, weight='bold', fontsize=30)
ax2.bar(x=x_axis, height=percentageoftopABC, edgecolor='black', color='#4fd4db')
add_labels(x_axis, percentageoftopABC, axis=ax2)
ax2.set_title('World Population (% of variants recognised by TopABC)', loc='left', y=1, weight='bold', fontsize=30)
plt.xticks(fontsize=9, fontweight='bold', rotation='vertical')
plt.yticks(fontsize=9, fontweight='bold')

plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
            'Recognition_peptides_TopABC_NCBI.pdf', dpi=450)
plt.show()

# %% ----------------------------------------------------------------

newdf = Dataframe.drop(['Number of genomes'], axis=1)
newdf = newdf.drop(['HLAs recognised'], axis=1)
# newdf = newdf.filter(items=list_topABC_2)
row_sums = newdf.apply(lambda x: x.sum())
newdf = newdf.transpose()
newdf['Sum'] = row_sums
fig, ax = plt.subplots(figsize=(15, 10))
sns.barplot(newdf, x=newdf.index, y='Sum')
add_labels(x_axis, percentageoftopABC, axis=ax2)
plt.xticks(fontsize=9, fontweight='bold', rotation='vertical')
plt.yticks(fontsize=9, fontweight='bold')
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
            'HLA_Peptide_recognition_all_lengths.pdf', dpi=600)
plt.show()

# %% ----------------------------------------------------------------

# Read the data from the Excel files.
length8 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                        'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
                        'TopABC_NCBI_length8_summarydata.xlsx',
                        sheet_name='summarydata', index_col=0)
length9 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                        'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
                        'TopABC_NCBI_length9_summarydata.xlsx',
                        sheet_name='summarydata', index_col=0)
length10 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                         'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
                         'TopABC_NCBI_length10_summarydata.xlsx',
                         sheet_name='summarydata', index_col=0)
length11 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                         'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
                         'TopABC_NCBI_length11_summarydata.xlsx',
                         sheet_name='summarydata', index_col=0)

length8 = length8.drop('Number of genomes', axis=1).drop(['HLAs recognised'], axis=1)
# length8 = length8.filter(items=list_topABC_2)
length9 = length9.drop('Number of genomes', axis=1).drop(['HLAs recognised'], axis=1)
# length9 = length9.filter(items=list_topABC_2)
length10 = length10.drop('Number of genomes', axis=1).drop(['HLAs recognised'], axis=1)
# length10 = length10.filter(items=list_topABC_2)
length11 = length11.drop('Number of genomes', axis=1).drop(['HLAs recognised'], axis=1)
# length11 = length11.filter(items=list_topABC_2)

# Transpose all datasets.
length8 = length8.T
length9 = length9.T
length10 = length10.T
length11 = length11.T

# Plot the data.
fig, axes = plt.subplots(4, 1, figsize=(25, 25), sharex=False)

sns.barplot(x=length8.index, y=length8.sum(axis=1), ax=axes[0])
axes[0].set_title('Length 8')

sns.barplot(x=length9.index, y=length9.sum(axis=1), ax=axes[1])
axes[1].set_title('Length 9')

sns.barplot(x=length10.index, y=length10.sum(axis=1), ax=axes[2])
axes[2].set_title('Length 10')

sns.barplot(x=length11.index, y=length11.sum(axis=1), ax=axes[3])
axes[3].set_title('Length 11')

for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/NCBI_CSP/resultsPredictionBinding_NCBI_only_TopABC/'
            'HLAs_peptide_recognition_lengths_divided.pdf', dpi=1000)
plt.show()

# %% New dictionary of HLAs per population

Table = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/FilteredDataAllele.xlsx',
                      sheet_name='TOP_ABC')

result_dict = Table.to_dict(orient='dict')

new_df = Dataframe.copy().drop('Number of genomes', axis=1).drop(['HLAs recognised'], axis=1)
data = {'x_axis': [], 'y_axis': []}
for country in result_dict.keys():
    alleles = [str(x) for x in result_dict[country].values()]
    # topABC_2 = set(list_topABC_2)
    # topABC_2_intersect = topABC_2.intersection(set(alleles))
    country_summary = new_df.filter(alleles)
    hits = np.sum(country_summary.values)
    data['x_axis'].append(country)
    data['y_axis'].append(hits)

fig, ax = plt.subplots(figsize=(25, 20))
# Creation of the plot with the different populations
sns.barplot(x='x_axis', y='y_axis', data=data, edgecolor='black')
ax.set_title('Number of peptides recognised by the Top ABC per country', loc='left', y=1, weight='bold', fontsize=30)
plt.xticks(fontsize=20, fontweight='bold', rotation='vertical')
plt.yticks(fontsize=20, fontweight='bold')

plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/NCBI_CSP/'
            'resultsPredictionBinding_NCBI_only_TopABC/Population_peptide_recognition.pdf', dpi=600)
plt.show()

# %% Preparing the different HLAs prevalence in different populations
