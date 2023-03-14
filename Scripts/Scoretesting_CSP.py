#!/usr/bin/env python3

# Install all packages before running main script, look Documentation file on how to do it.
import os

import pandas as pd
from Bio import SeqIO

if __name__ == '__main__':
    """ 
    It is important to have prepared the files in a folder first, the following script uses the table produced by
    the software that tests the binding of MHC-class I.

    This script is run by typing: 

    """

    print('Please introduce the name of the folder where you have store the data to analyse:')

    wd = input()

    os.chdir(str(wd))  # Change working directory

    # Import arguments into variables.
    print('Please introduce the name of the file where you have store the data to analyse')

    inputfile = input()
    MHCtable = pd.read_excel(inputfile, sheet_name='Sheet1')

    print('Please introduce the name of the fasta file')

    inputfasta = input()
    fastaprotein = SeqIO.parse(str(inputfasta), 'fasta')

    print('Please introduce the threshold')
    number = input()

    MHCtable = MHCtable[MHCtable['percentile_rank'] < int(number)]

    Good_Alleles = ['HLA-A*01:01', 'HLA-B*08:01']

    for record in fastaprotein:
        sequence = str(record.seq)
        for begining_pep in range(len(sequence)):
            end_pep = begining_pep + 10
            if end_pep <= len(sequence):
                counter = 0

                for index, row in MHCtable.iterrows():
                    start = int(MHCtable.loc[index]['start'])
                    end = int(MHCtable.loc[index]['end'])
                    if (start >= begining_pep) and (end <= end_pep):
                        counter = counter + 1
                    else:
                        continue

                print(
                    'Peptide between ' + str(begining_pep) + ' and ' + str(end_pep) + ' contains ' + str(counter) + ' hits')
                with open('Output.tsv', 'a') as f:
                    f.write('Peptide' + "\t" + str(begining_pep) + "\t" + str(end_pep) + "\t" + str(counter) + "\n")

