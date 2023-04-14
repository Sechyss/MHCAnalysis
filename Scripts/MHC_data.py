import argparse
import textwrap

import pandas as pd


# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Analysis of Matlab simulations results
------------------------------------------
Tasks added so far:
allHLAs       -- Output all the HLAs for that particular population group
filteredHLAs  -- Output filtered HLAs for that particular population group taking into account those already added

'''))


parser.add_argument("-t", "--task", dest="task", help="Task to perform", type=str, choices=["allHLAs",
                                                                                            "filteredHLAs",
                                                                                            'merge-tables'])
parser.add_argument("-f", "--freqtable", metavar='table.xlsx', dest="allelefreqtable", help="Allele frequency table",
                    type=str)
parser.add_argument("-l", "--list", metavar='list.txt', dest="listhla", help="List of HLA in algorithm", type=str)
parser.add_argument("-p", "--population", dest="population", help="String Population (all pops)", type=str)
parser.add_argument("-x", "--spreadsheet", dest="spreadsheet", help="Excel spreadsheet with data", type=str)
parser.add_argument("-y", "--spreadsheet", dest="spreadsheet2", help="Excel spreadsheet with data", type=str)
parser.add_argument("-o", "--output", dest="output", help="Output file", type=str)

args = parser.parse_args()


def all_hlas():
    df = pd.read_excel(args.allelefreqtable, sheet_name='Sheet1', index_col=0)
    relative_allele_freq = df[df['Population'].str.contains(args.population)]

    list_alleles = set(relative_allele_freq['Allele'])

    mhc_predict_df = pd.read_table(args.listhla, header=None, sep='\t')

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

    print(allele)
    return finaldf


def filtered_hla():
    alleles_df = all_hlas()
    excel = pd.read_excel(args.spreadsheet, sheet_name='summarydata', index_col=0)
    list_to_filter = list(excel.columns)
    list_to_filter.pop(0)
    filtered_df = alleles_df[~alleles_df[1].isin(list_to_filter)]
    allele = list(filtered_df[1])
    length = list(filtered_df[2])

    print(','.join(allele), ','.join(map(str, length)))


def combinetables():
    excel = pd.read_excel(args.spreadsheet, sheet_name='summarydata', index_col=0)
    excel2 = pd.read_excel(args.spreadsheet2, sheet_name='summarydata', index_col=0)
    excel3 = pd.merge(left=excel, right=excel2, how='outer', left_index=True, right_index=True)
    excel3 = excel3.reindex(sorted(excel3.columns), axis=1)
    excel3.to_excel(args.output, sheet_name='summarydata')


def main():
    task = {'allHLAs': all_hlas, 'filteredHLAs': filtered_hla, 'merge-tables': combinetables}
    choosetask = str(args.task)
    task[choosetask]()


if __name__ == '__main__':
    main()
