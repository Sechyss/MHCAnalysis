import argparse
import textwrap

import pandas as pd

# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Printing of HLAs to run the prediction tool.
------------------------------------------
Tasks added so far:
allHLAs       -- Output all the HLAs for that particular population group
filteredHLAs  -- Output filtered HLAs for that particular population group taking into account those already added
merge-tables  -- Merge tables from different results into one single spreadsheet

'''))

# Parse the input for the different scripts
parser.add_argument("-f", "--freqtable", metavar='table.xlsx', dest="allelefreqtable", help="Allele frequency table",
                    type=str)
parser.add_argument("-l", "--coordinates", metavar='coordinates.txt', dest="listhla", help="List of HLA in algorithm", type=str)
parser.add_argument("-p", "--population", dest="population", help="String Population (all pops)", type=str)
parser.add_argument("-x", "--spreadsheet", dest="spreadsheet", help="Excel spreadsheet with data", type=str)
parser.add_argument("-y", "--spreadsheet2", dest="spreadsheet2", help="Excel spreadsheet with data", type=str)
parser.add_argument("-o", "--output", dest="output", help="Output file", type=str)
parser.add_argument("-t", "--task", dest="task", help="Task to perform", type=str, choices=["allHLAs",
                                                                                            "filteredHLAs",
                                                                                            "merge-tables"])

# Compile the coordinates of flags for the script
args = parser.parse_args()


def all_hlas():
    df = pd.read_excel(args.allelefreqtable, sheet_name='Combination_HLAs', index_col=0)
    relative_allele_freq = df.filter(like=args.population)
    # set the threshold
    threshold = 0.01

    # filter the DataFrame to select only rows whose values are over the threshold for all columns
    filtered_df = relative_allele_freq[relative_allele_freq.apply(
        lambda row: all(val >= threshold for val in row), axis=1)]

    list_alleles = set(filtered_df.index)

    mhc_predict_df = pd.read_table(args.listhla, header=None, sep='\t')

    # Replace the empty spaces in the coordinates
    mhc_predict_df[1] = mhc_predict_df[1].apply(lambda x: str(x).replace(' ', ''))

    # Replace part of the strings to see if they match with the other lists
    list_smm = list(mhc_predict_df[1])
    temp_smm = [x.replace('HLA-', '').replace(' ', '') for x in list_smm]
    common_mhc = list(list_alleles & set(temp_smm))

    # Regain the HLA- part of the string to print correct output
    finalist_common = ['HLA-' + str(x) for x in common_mhc]

    finaldf = mhc_predict_df[mhc_predict_df[1].isin(finalist_common)]

    # This next part of the code can be tweaked to select a specific kmer size

    finaldf = finaldf[finaldf[2] <= 11]  # <-- Change this line to select a specific kmer size

    allele = list(finaldf[1])
    length = list(finaldf[2])
    # Print the coordinates of alleles and lengths to copy into the command line for prediction tool
    print(','.join(allele), ','.join(map(str, length)))

    return finaldf


def filtered_hla():
    # Extract all the alleles that match the coordinates using the previous function all_hlas()
    alleles_df = all_hlas()
    # Extract the alleles from the table containing previously studied alleles
    excel = pd.read_excel(args.spreadsheet, sheet_name='summarydata', index_col=0)
    list_to_filter = list(excel.columns)
    list_to_filter.pop(0)

    # Filter those already included to avoid duplicates and wasting time and memory
    filtered_df = alleles_df[~alleles_df[1].isin(list_to_filter)]
    allele = list(filtered_df[1])
    length = list(filtered_df[2])

    # Print the coordinates of alleles that haven't been filtered for prediction tool
    print(','.join(allele), ','.join(map(str, length)))

    return filtered_df


def combinetables():
    excel = pd.read_excel(args.spreadsheet, sheet_name='summarydata', index_col=0)
    excel2 = pd.read_excel(args.spreadsheet2, sheet_name='summarydata', index_col=0)
    excel3 = pd.merge(left=excel, right=excel2, how='outer', left_index=True, right_index=True)
    excel3 = excel3.reindex(sorted(excel3.columns), axis=1)
    excel3.to_excel(args.output, sheet_name='summarydata')


def main():
    task = {'allHLAs': all_hlas, 'filteredHLAs': filtered_hla, 'merge-tables': combinetables}

    # Read the flag argument for the task and call the function into main
    choosetask = str(args.task)
    task[choosetask]()


if __name__ == '__main__':
    main()
