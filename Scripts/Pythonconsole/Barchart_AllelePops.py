import pandas as pd
from matplotlib import pyplot as plt
import pickle
import seaborn as sns


def extractionofpercentages(dictionary):
    """
      Function will calculate the number of HLAs that recognise the Kmers

      Args:
        dictionary: The dictionary to slice.

      Returns:
        A new list of percentages to plot later
      """
    #
    list_percentages = []  # Empty object to store the information
    values = list(dictionary.values())  # Extraction of the values from the dictionary
    total_values = list(set(item for sublist in values for item in sublist))  # Flatten the list of values
    total_values.sort()  # Sort the list in order
    counter = 0
    for ids in number_variants:
        list_ids = range(counter, ids + counter)  # Create the list of ranges which represent the IDs
        hits = 0
        for seqid in list_ids:  # Run through the IDs
            if seqid in total_values:  # If an ID is found in an HLA list it will add up to the hit counter
                hits += 1
        percentage_of_hits = round((hits / ids * 100), 2)  # Normalize based on the total number of possibilities
        list_percentages.append(percentage_of_hits)
        counter += ids  # This part will limit future ranges of IDs

    return list_percentages


def add_labels(x, z, axis):  # Add labels to the plot later
    for i in range(len(x)):
        axis.text(i, z[i], z[i], ha='center', fontweight='bold', color='black', fontsize=5)


def slice_dict(dictionary, keys):
    """
  Slices a dictionary to extract only keys and values that match a list of keys.

  Args:
    dictionary: The dictionary to slice.
    keys: The list of keys to match.

  Returns:
    A new dictionary containing only the keys and values that match the list of keys.
  """
    new_dict = {}
    for key in keys:
        if key in dictionary:
            new_dict[key] = dictionary[key]
    return new_dict


# %% Upload the data and the functions to run the plot
Dataframe = pd.read_excel('/Users/u2176312/OneDrive - University of '
                          'Warwick/CSP/AllelePops/Kmer_CSP_region_273-375_TopABC_summarydata_data.xlsx',
                          sheet_name='summarydata', index_col=0)
tempfile = open('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                '/AllelePops/Allele_seq_ids_273-375_region_RTSS_alllengths.pickle', 'rb')
dictionary_alleles = pickle.load(tempfile)
x_axis = Dataframe.index.to_numpy()
number_variants = Dataframe['Variants'].to_numpy()

# %%
hits = extractionofpercentages(dictionary_alleles)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
# Creation of the plot with the different populations
ax1.bar(x=x_axis, height=number_variants, edgecolor='black', color='white')
add_labels(x_axis, number_variants, axis=ax1)
ax1.set_title('Number of variants', loc='left', y=1, weight='bold', fontsize=30)
ax2.bar(x=x_axis, height=hits, edgecolor='black', color='#4fd4db')
add_labels(x_axis, hits, axis=ax2)
ax2.set_title('World Population (% of variants recognised)', loc='left', y=1, weight='bold', fontsize=30)
plt.xticks(fontsize=9, fontweight='bold', rotation='vertical')
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/AllelePops/TopABCResults/Kmer_Variant_recognition_topABC.pdf', dpi=300)
plt.show()

# %% ----------------------------------------------------------------

newdf = Dataframe.drop(['Variants'], axis=1)
row_sums = newdf.apply(lambda x: x.sum())
newdf = newdf.transpose()
newdf['Sum'] = row_sums
fig, ax = plt.subplots(figsize=(15, 10))
sns.barplot(newdf, x=newdf.index, y='Sum')
add_labels(x_axis, hits, axis=ax2)
plt.xticks(fontsize=9, fontweight='bold', rotation='vertical')
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/AllelePops/TopABCResults/HLA_Variant_recognition_alllenghts.pdf', dpi=300)
plt.show()

# %% ----------------------------------------------------------------

# Read the data from the Excel files.
length8 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                        'Warwick/CSP/AllelePops/Kmer_CSP_region_273-375_summarydata_length8.xlsx',
                        sheet_name='summarydata', index_col=0)
length9 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                        'Warwick/CSP/AllelePops/Kmer_CSP_region_273-375_summarydata_length9.xlsx',
                        sheet_name='summarydata', index_col=0)
length10 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                         'Warwick/CSP/AllelePops/Kmer_CSP_region_273-375_summarydata_length10.xlsx',
                         sheet_name='summarydata', index_col=0)
length11 = pd.read_excel('/Users/u2176312/OneDrive - University of '
                         'Warwick/CSP/AllelePops/Kmer_CSP_region_273-375_summarydata_length11.xlsx',
                         sheet_name='summarydata', index_col=0)

length8 = length8.drop('Variants', axis=1)
length9 = length9.drop('Variants', axis=1)
length10 = length10.drop('Variants', axis=1)
length11 = length11.drop('Variants', axis=1)

# Transpose all datasets.
length8 = length8.T
length9 = length9.T
length10 = length10.T
length11 = length11.T

# Plot the data.
fig, axes = plt.subplots(4, 1, figsize=(20, 10), sharex=False)

sns.barplot(x=length8.index, y=length8.sum(axis=1), ax=axes[0])
axes[0].set_title('Length 8')

sns.barplot(x=length9.index, y=length9.sum(axis=1), ax=axes[1])
ax2.set_title('Length 9')

sns.barplot(x=length10.index, y=length10.sum(axis=1), ax=axes[2])
axes[2].set_title('Length 10')

sns.barplot(x=length11.index, y=length11.sum(axis=1), ax=axes[3])
axes[3].set_title('Length 11')

for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/AllelePops/TopABCResults/HLAs_Variant_recognition_lengths_divided.pdf', dpi=300)
plt.show()

# %%% Plot divided by countries

temp_file = open('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Dictionary_country_TopABC.pickle',
                 'rb')
dict_alleles = pickle.load(temp_file)
list_countries = [x for x in dict_alleles.keys()]
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
new_df['Variants'] = Dataframe['Variants']

for country in dict_alleles.keys():
    alleles = [str(x) for x in dict_alleles[country]]
    country_summary = slice_dict(dictionary_alleles, alleles)
    hits = extractionofpercentages(country_summary)
    new_df[country] = hits
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(25, 20), sharex=True)
    # Creation of the plot with the different populations
    ax1.bar(x=x_axis, height=number_variants, edgecolor='black', color='white')
    ax1.set_title('Number of variants', loc='left', y=1, weight='bold', fontsize=30)
    ax2.bar(x=x_axis, height=hits, edgecolor='black', color=color_dictionary[country])
    ax2.set_title(country + ' (% of variants recognised)', loc='left', y=1, weight='bold', fontsize=30)

    plt.xticks(fontsize=10, fontweight='bold', rotation='vertical')
    # Get the yticks for the subplots
    #yticks1 = ax1.get_yticks()
    #yticks2 = ax2.get_yticks()

    # Increase the fontsize and fontweight of the yticks
    #yticks1 = [x.set_fontsize(15).set_fontweight('bold') for x in yticks1]
    #yticks2 = [x.set_fontsize(15).set_fontweight('bold') for x in yticks2]
    # Set the yticks for the subplots
    #ax1.set_yticks(yticks1)
    #ax2.set_yticks(yticks2)
    plt.tight_layout()
    plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/TopABCResults/CountryFigures/'
                'Kmer_Variant_recognition_' + country + '.pdf', dpi=300)
    plt.show()

new_df.to_csv(
    '/Users/u2176312/OneDrive - University of '
    'Warwick/CSP/AllelePops/TopABCResults/Kmer_Variant_recognition_world_rawdata.csv')

