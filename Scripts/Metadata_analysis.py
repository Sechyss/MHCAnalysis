# Upload of packages to run the script
import argparse
import textwrap
import warnings

import pandas as pd
import pysam
import seaborn as sns
from matplotlib import pyplot as plt

warnings.simplefilter(action='ignore', category=FutureWarning)
# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Analysis of metadata from VCF file and all of its variants.
-------------------------------------------------------------
'''))
# Parse command line arguments
# -i FILE -t TASK -h HELP -o OUTPUT -th Threshold

parser.add_argument("-v", "--vcf", metavar='file.vcf.gz', dest="vcf", help="VCF file", type=str)
parser.add_argument("-m", "--mhc", metavar='mhc', dest="mhc", type=str)
parser.add_argument("-t", "--table", metavar='file.xlsx', dest="table", help="Excel spreadsheet", type=str)

args = parser.parse_args()

# =============================================================================

metadata = pd.read_excel(args.table, sheet_name='Sheet1', index_col=0)
vcf_reader = pysam.VariantFile(str(args.vcf), 'r')
mhc_list = pd.read_excel(args.mhc, sheet_name='Sheet1', index_col=0)
higher_freq = mhc_list[mhc_list['Allele_freq'] > 0.01]
countries = higher_freq['Region'].unique().tolist()

fig, ax = plt.subplots()

for country in countries:
    sns.color_palette("Paired", len(countries))
    data = higher_freq[higher_freq['Region'] == country]
    sns.distplot(data['Allele_freq'], kde='True', label=country)

plt.legend(loc='best', frameon=False, prop={'size': 7})
plt.show()

for record in vcf_reader.fetch():
    samples = record.samples.keys()
    break

# "Venn diagram" of MHCs in each country

collecting_dict = {}
for country in countries:
    mhcs_per_country_df = mhc_list[mhc_list['Region'] == country]
    list_to_append = list(set(mhcs_per_country_df['Allele']))
    list_to_append.sort()
    collecting_dict.update({country: list_to_append})

df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in collecting_dict.items()]))
