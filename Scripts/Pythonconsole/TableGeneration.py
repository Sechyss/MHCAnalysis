import requests
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np


def extract_table_from_html(url):
    html = requests.get(url).content
    df_list = pd.read_html(html)
    df = df_list[-1]  # the last table in the list is usually the one we want
    df.drop(['Line', '% of individuals that have the allele'], axis=1, inplace=True)
    return df


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

new_data = pd.DataFrame()

#writer = pd.ExcelWriter('/Users/u2176312/OneDrive - University of Warwick/CSP/RawDataAllele_freq.xlsx',
#                        engine='openpyxl')
for i in tqdm(range(len(list_populations))):
    final_df = extract_table_from_html(urls[i])
    new = final_df.filter(['Allele', 'Allele Frequency'], axis=1)
    new.set_index('Allele', inplace=True)
    new = new.rename(columns={'Allele Frequency': str(list_populations[i])})
    new_data = pd.merge(new_data, new, how='outer', left_index=True, right_index=True)
#    final_df.to_excel(writer, sheet_name=list_populations[i])
new_data = new_data.drop(index=new_data[new_data.index.str.contains('D')].index)
new_data = new_data.fillna(0)
#new_data.to_excel(writer, sheet_name='Combination_HLAs')

#writer.close()

# %% Exploration of the data to study a threshold

new_data['mean'] = new_data.mean(axis=1)
plt.figure(figsize=(10, 5))
plt.xlabel('Allele frequency')
plt.ylabel('Counts')
counts, bins = np.histogram(new_data['mean'])
plt.stairs(counts, bins)
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')
plt.tight_layout()
plt.show()
