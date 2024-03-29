import requests
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import numpy as np
import itertools
import pickle

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def extract_table_from_html(url):  # Function to extract table from html address
    html = requests.get(url).content
    df_list = pd.read_html(html)
    df = df_list[-1]  # the last table in the coordinates is usually the one we want
    df.drop(['Line', '% of individuals that have the allele'], axis=1, inplace=True)
    return df

# Creation of the coordinates of populations and their url to download. Important to check the order of the lists
list_populations = ['Papua New Guinea Karimui Plateau Pawaia',
                    'Papua New Guinea Madang',
                    'New Caledonia',
                    'Philippines Ivatan',
                    'USA Hawaii Okinawa',
                    'China Guizhou Province Shui',
                    'India Khandesh Region Pawra',
                    'Georgia Tibilisi',
                    'Iran Gorgan',
                    'Mali Bandiagara',
                    'Morocco Nador Metalsa pop 2',
                    'Ghana Ga-Adangbe',
                    'Kenya Nyanza Province Luo tribe',
                    'Kenya Nandi',
                    'Zimbabwe Harare Shona',
                    'England North West',
                    'Finland',
                    'Greece pop 6',
                    'Ireland South'
                    ]
urls = [
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1737&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1714&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1712&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1231&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2058&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2408&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1297&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2026&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=3662&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1486&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1492&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=3108&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=3393&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=1484&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2057&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2837&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2028&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=3076&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard=',
    'http://allelefrequencies.net/hla6006a_scr.asp?hla_locus=&hla_locus_type=Classical&hla_allele1=&hla_allele2=&hla_selection=&hla_pop_selection=&hla_population=2275&hla_country=&hla_dataset=&hla_region=&hla_ethnic=&hla_study=&hla_sample_size=&hla_sample_size_pattern=&hla_sample_year=&hla_sample_year_pattern=&hla_level=&hla_level_pattern=&hla_show=&hla_order=order_1&standard='
]

# Creation of the datframe and dictionary for further analysis. Dataframe will be saved in a Excel file
new_data = pd.DataFrame()

#writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/RawDataAllele_freq.xlsx',
#                        engine='openpyxl')
dict_alleles = {}
for i in tqdm(range(len(list_populations))):
    final_df = extract_table_from_html(urls[i])
    new = final_df.filter(['Allele', 'Allele Frequency'], axis=1)
    listalleles= final_df['Allele'].tolist()
    filteredMHCI = [x for x in listalleles if 'D' not in x]
    dict_alleles.update({str(list_populations[i]): filteredMHCI})
    new.set_index('Allele', inplace=True)
    new = new.rename(columns={'Allele Frequency': str(list_populations[i])})
    new_data = pd.merge(new_data, new, how='outer', left_index=True, right_index=True)
#    final_df.to_excel(writer, sheet_name=list_populations[i])
new_data = new_data.drop(index=new_data[new_data.index.str.contains('D')].index)
new_data = new_data[new_data.index.str.count(':') < 2]
new_data = new_data.fillna(0)
#new_data.to_excel(writer, sheet_name='Combination_HLAs')

#writer.close()
with open('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Dictionary_country_alleles.pickle', 'wb') as f:
    pickle.dump(dict_alleles, f)
# %% Exploration of the data to study a threshold

new_data = new_data[new_data.max(axis=1) > 0.15]
descriptivedata = new_data.transpose().replace(0, np.nan).describe()

#%% Plotting the data and exploration
#descriptivedata.to_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/AlleleFreqdescriptive.csv')

plt.figure(figsize=(10,6))
plt.xlabel('Allele frequency', fontweight='bold', fontsize=15)
plt.ylabel('Counts', fontweight='bold', fontsize=15)
sns.histplot(data=new_data.replace(0, np.nan), legend=True, bins=40)
plt.xticks(fontweight='bold')
plt.locator_params(axis='x', nbins=30)
plt.yticks(fontweight='bold')
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Allele_frequencies_histogram.pdf', dpi=300)
plt.show()


plt.figure(figsize=(95, 30))
sns.boxplot(data=new_data.transpose().replace(0, np.nan))
plt.xlabel('Alleles', fontsize=75, fontweight='bold')
plt.ylabel('Allele frequency', fontsize=75, fontweight='bold')
plt.yticks(fontweight='bold', fontsize=50)
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Allele_frequencies_boxplot_allele.pdf', dpi=300)
plt.show()

plt.figure(figsize=(10, 10))
sns.boxplot(data=new_data.replace(0, np.nan))
plt.xlabel('Populations', fontweight='bold', fontsize=15)
plt.ylabel('Allele frequency', fontweight='bold', fontsize=15)
plt.locator_params(axis='y', nbins=34)
plt.yticks(fontweight='bold')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Allele_frequencies_boxplot_pops.pdf', dpi=300)
plt.show()


filteredindex = new_data.index.tolist()
x_values = dict_alleles.keys()
y_values = []
for x in x_values:
    listvalues = dict_alleles[x]
    inter = set(listvalues) & set(filteredindex)
    y_values.append(len(inter))

data = pd.DataFrame.from_dict(dict(zip(x_values, y_values)), orient='index', columns=['Number of Alleles'])


plt.figure(figsize=(10, 10))
sns.barplot(data=data.transpose(), edgecolor='black')
plt.xlabel('Populations', fontweight='bold', fontsize=15)
plt.ylabel('Number of HLAs I', fontweight='bold', fontsize=15)
plt.yticks(fontweight='bold')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Allele_number_pops.pdf', dpi=300)
plt.show()


shared = {}


for key1 in dict_alleles:
    for key2 in dict_alleles:
        if key1 != key2:
            alleles1 = [x for x in dict_alleles[key1] if x in filteredindex]
            alleles2 = [x for x in dict_alleles[key2] if x in filteredindex]
            intersection = set(alleles1) & set(alleles2)
            if intersection:
                shared[(key1, key2)] = len(intersection) / len(set(alleles1 + alleles2))

keys = sorted(list(set([key[0] for key in shared.keys()] + [key[1] for key in shared.keys()])))
data1 = [[shared.get((key1, key2), 0) for key2 in keys] for key1 in keys]

# Create a mask for the upper triangle of the heatmap
mask = np.triu(np.ones_like(data1, dtype=bool))

plt.figure(figsize=(20, 15))
sns.heatmap(data1, xticklabels=keys, yticklabels=keys, annot=True,
            cbar=False, cmap="YlGnBu", mask=mask, vmin=0, vmax=1, fmt=".2f")
plt.yticks(fontweight="bold")
plt.xticks(fontweight="bold")

# add colorbar legend
norm = colors.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap="YlGnBu", norm=norm)
sm.set_array([])

cbar = plt.colorbar(sm)
cbar.ax.set_ylabel('Percentage of Shared Values', rotation=270, labelpad=15)

plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Allele_sharing_heatmap.pdf', dpi=300)
plt.show()

#%% Print the table and coordinates of alleles and their lenghths
new_data = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/FilteredDataAllele.xlsx',
                         sheet_name='FilteredAlleles', index_col=0)

# new_data = new_data[new_data.max(axis=1) > 0.01]

listallelesnetpanmhc = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/ListMHC_humans_netmhcpan.txt',
                               header=None, sep='\t')


All_HLAs = new_data.index.tolist()
temp_hlas = ['HLA-' + str(x) for x in All_HLAs]
netmhcalleles = list(listallelesnetpanmhc[0])
temp_list = [x.replace(' ', '') for x in netmhcalleles]
common_mhc = list(set(temp_hlas) & set(temp_list))
print('" "'.join(set(common_mhc)))

lengths = [8,9,10,11]

combinations = list(itertools.product(All_HLAs, lengths))

alleles = ['HLA-' + str(x[0]) for x in combinations]
lengths = [x[1] for x in combinations]

#print(','.join(alleles), ','.join(map(str, lengths)))


# %% Sanity Checks fpr the full coordinates vs successful ones

listsucdess_all = pd.read_excel('/Users/u2176312/OneDrive - University of '
                          'Warwick/CSP/AllelePops/All_HLAs_Kmer_CSP_region_273-375_summarydata_length_all.xlsx',
                          sheet_name='summarydata', index_col=0)
list_success = listsucdess_all.columns.tolist()
list_success.remove('Variants')

intersection_vals = list(set(common_mhc).difference(set(list_success)))
print('" "'.join(set(intersection_vals[:])))

#%% Print top ABCs

new_data = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/FilteredDataAllele.xlsx',
                         sheet_name='TOP_ABC')
values = new_data.values.ravel()
flattenlist = set(list(values))
listallelesnetpanmhc = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/ListMHC_humans_netmhcpan.txt',
                               header=None, sep='\t')
netmhcalleles = list(listallelesnetpanmhc[0])
temp_list = [x.replace(' ', '') for x in netmhcalleles]
common_mhc = list(flattenlist & set(temp_list))
print('" "'.join(set(common_mhc)))

new_dict = {}
for column in new_data.columns:
    new_dict.update({column: new_data[column].tolist()})

with open('/Users/u2176312/OneDrive - University of Warwick/CSP/AllelePops/Dictionary_country_TopABC.pickle', 'wb') as f:
    pickle.dump(new_dict, f)