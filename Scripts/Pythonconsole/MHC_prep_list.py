#!/usr/bin/env/MHC python3
import pandas as pd

df = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/SubSaharaPop_allelefreq.xlsx',
                   sheet_name='Sheet1', index_col=0)
relative_allefreq = df[df['Population'] == 'Kenya']

list_alleles = set(relative_allefreq['Allele'])

mhc_predict_df = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/ListMHC_humans_smm.txt',
                               header=None, sep='\t')

mhc_predict_df[1] = mhc_predict_df[1].apply(lambda x: str(x).replace(' ', ''))

list_smm = list(mhc_predict_df[1])
temp_smm = [x.replace('HLA-', '').replace(' ', '') for x in list_smm]
common_mhc = list(list_alleles & set(temp_smm))

finalist_common = ['HLA-' + str(x) for x in common_mhc]

finaldf = mhc_predict_df[mhc_predict_df[1].isin(finalist_common)]

# This next part of the code can be removed to select a specific kmer size

finaldf = finaldf[finaldf[2] <= 12]  # <-- Change this line to select a specific kmer size


allele = list(finaldf[1])
length = list(finaldf[2])
print(','.join(allele), ','.join(map(str, length)))
