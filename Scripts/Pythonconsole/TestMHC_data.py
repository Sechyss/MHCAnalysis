import pandas as pd

df = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/SubSaharaPop_allelefreq.xlsx',
                   sheet_name='Sheet1', index_col=0)
relative_allefreq = df[df['Population'].str.contains('Kenya')]

list_alleles = set(relative_allefreq['Allele'])

mhc_predict_df = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/ListHLA_netmhcpan_el.txt',
                               header=None, sep='\t')

# Replace the empty spaces in the list
mhc_predict_df[1] = mhc_predict_df[1].apply(lambda x: str(x).replace(' ', ''))

# Replace part of the strings to see if they match with the other lists
list_smm = list(mhc_predict_df[1])
temp_smm = [x.replace('HLA-', '').replace(' ', '') for x in list_smm]
common_mhc = list(list_alleles & set(temp_smm))

# Regain the HLA- part of the string to print correct output
finalist_common = ['HLA-' + str(x) for x in common_mhc]

finaldf = mhc_predict_df[mhc_predict_df[1].isin(finalist_common)]

# This next part of the code can be removed to select a specific kmer size

finaldf = finaldf[finaldf[2] <= 11]  # <-- Change this line to select a specific kmer size

allele = list(finaldf[1])
length = list(finaldf[2])
print(','.join(allele), ','.join(map(str, length)))

excel = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/CSP_SNP_region/'
                      'Kmer_CSP_region_283-302_summarydata_Finland.xlsx', sheet_name='summarydata', index_col=0)
filtered_hla = list(excel.columns)
filtered_hla.pop(0)

filtered_df = finaldf[~finaldf[1].isin(filtered_hla)]
notfiltered = finaldf[finaldf[1].isin(filtered_hla)]
allele = list(filtered_df[1])
length = list(filtered_df[2])
allele2 = list(notfiltered[1])
length2 = list(notfiltered[2])

print(','.join(allele), ','.join(map(str, length)))
print(', '.join(allele2), ','.join(map(str, length2)))

excel2 = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/CSP_SNP_region/'
                       'Kmer_CSP_region_283-307_summarydata_Kenya.xlsx', sheet_name='summarydata', index_col=0)
excel3 = pd.concat([excel, excel2], join='outer', axis=1)


def same_merge(x): return ','.join(x[x.notnull()].astype(str))


df_new = excel3.groupby(level=0, axis=1).apply(lambda x: x.apply(same_merge, axis=1))
