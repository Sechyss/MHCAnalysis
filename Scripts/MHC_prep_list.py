#!/usr/bin/env/MHC python3
import argparse
import textwrap

import pandas as pd

# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Summary table of the prediction binding protein in combination with the Sequences generated by the combination of all 
the SNPs in fasta file.
------------------------------------------
'''))

parser.add_argument("-f", "--freqtable", metavar='table.xlsx', dest="allelefreqtable", help="Allele frequency table",
                    type=str)
parser.add_argument("-l", "--list", metavar='list.txt', dest="listhla", help="List of HLA in algorithm", type=str)
parser.add_argument("-p", "--population", dest="population", help="String Population (all pops)", type=str)


args = parser.parse_args()

# =============================================================================

df = pd.read_excel(args.allelefreqtable, sheet_name='Sheet1', index_col=0)
relative_allefreq = df[df['Population'].str.contains(args.population)]

list_alleles = set(relative_allefreq['Allele'])

mhc_predict_df = pd.read_table(args.listhla, header=None, sep='\t')

mhc_predict_df[1] = mhc_predict_df[1].apply(lambda x: str(x).replace(' ', ''))

list_smm = list(mhc_predict_df[1])
temp_smm = [x.replace('HLA-', '').replace(' ', '') for x in list_smm]
common_mhc = list(list_alleles & set(temp_smm))

finalist_common = ['HLA-' + str(x) for x in common_mhc]

finaldf = mhc_predict_df[mhc_predict_df[1].isin(finalist_common)]

# This next part of the code can be removed to select a specific kmer size

finaldf = finaldf[finaldf[2] <= 11]  # <-- Change this line to select a specific kmer size


allele = list(finaldf[1])
length = list(finaldf[2])
print(','.join(allele), ','.join(map(str, length)))

print(allele)
