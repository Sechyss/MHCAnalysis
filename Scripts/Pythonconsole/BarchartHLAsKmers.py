import pandas as pd
from matplotlib import pyplot as plt
import pickle

Finland_region = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/CSP_SNP_region/'
                               'Kmer_CSP_region_283-302_summarydata_Finland.xlsx', sheet_name='summarydata',
                               index_col=0)
Kenya_region = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/CSP/CSP_SNP_region/'
                             'Kmer_CSP_region_283-307_summarydata_Kenya.xlsx', sheet_name='summarydata',
                             index_col=0)
tempfile = open('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                '/CSP_SNP_region/Allele_seq_ids.pickle', 'rb')
dictionary_alleles_kenya = pickle.load(tempfile)
tempfile.close()
tempfile = open('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                'CSP_SNP_region/Allele_seq_ids_283-302_region_Finland.pickle', 'rb')
dictionary_alleles_finland = pickle.load(tempfile)
tempfile.close()

Kenya_region = Kenya_region.filter(items=list(Finland_region.index), axis=0)

x_axis = Kenya_region.index.to_numpy()
number_variants = Kenya_region['Variants'].to_numpy()


def extractionofpercentages(dictionary):
    list_percentages = []
    counter = 0
    values = list(dictionary.values())
    total_values = list(set(item for sublist in values for item in sublist))
    for ids in number_variants:
        list_ids = range(counter + ids)
        totalids = len(list_ids)
        hits = 0
        for seqid in list_ids:
            if seqid in total_values:
                hits += 1
        percentage_of_hits = round((hits/totalids * 100), 2)
        list_percentages.append(percentage_of_hits)
        counter += ids
    return list_percentages


def add_labels(x, z, axis):
    for i in range(len(x)):
        axis.text(i, z[i], z[i], ha='center', fontweight='bold', color='black', fontsize=10)


hitsKenya = extractionofpercentages(dictionary_alleles_kenya)
hitsFinland = extractionofpercentages(dictionary_alleles_finland)

# noinspection PyTypeChecker
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 12), sharex=True)

ax1.bar(x=x_axis, height=number_variants, edgecolor='black', color='white')
add_labels(x_axis, number_variants, axis=ax1)
ax1.set_title('Number of variants', loc='left', y=1, weight='bold')
ax2.bar(x=x_axis, height=hitsKenya, edgecolor='black', color='green')
add_labels(x_axis, hitsKenya, axis=ax2)
ax2.set_title('Kenya Population', loc='left', y=1, weight='bold')
ax3.bar(x=x_axis, height=hitsFinland, edgecolor='black', color='red')
add_labels(x_axis, hitsFinland, axis=ax3)
ax3.set_title('Finland Population', loc='left', y=1, weight='bold')

plt.xticks(fontsize=9, fontweight='bold')
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/CSP/CSP_SNP_region/'
            'ComparisonRegion283-302_Kenya_Finland_Population.pdf', dpi=300)

plt.show()
