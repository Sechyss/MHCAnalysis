#!/usr/bin/env python3

# Install all packages before running main script, look Documentation file on how to do it.
import os
import warnings

import matplotlib
import matplotlib.patches as mpatches
import pandas as pd
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord

# Pandas are deprecating append method but for the moment is the fastest way to append rows to a dataframe
warnings.simplefilter(action='ignore', category=FutureWarning)
# Setting the font for plots to be Times New Roman (can be changed)
matplotlib.rcParams['font.family'] = "Times New Roman"

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

    Good_Alleles = ['HLA-A*01:01', 'HLA-B*08:01']

    for record in fastaprotein:
        features = [GraphicFeature(start=1, end=18, color="olivedrab", label='Signal'),
                    GraphicFeature(start=93, end=97, color="brown", label='Region I'),
                    GraphicFeature(start=105, end=272, color="skyblue", label='Tandem repeats'),
                    GraphicFeature(start=322, end=375, color="crimson", label='TSP type-1'),
                    GraphicFeature(start=375, end=397, color="purple", label='GPI-anchor')]
        for index, row in MHCtable.iterrows():
            if int(MHCtable.loc[index]['percentile_rank']) <= int(number):
                allele = str(MHCtable.loc[index]['variant'])
                start = int(MHCtable.loc[index]['start'])
                end = int(MHCtable.loc[index]['end'])

                # Comparison with previous literature and collection of features for the plots
                if allele in Good_Alleles:
                    features.append(GraphicFeature(start=start, end=end, color="#ffd700",
                                                   label=allele))
                else:
                    features.append(GraphicFeature(start=start, end=end, color="#ccccff",
                                                   label=allele))

        # Plotting of the features on the protein sequence

        record = GraphicRecord(sequence_length=len(record.seq), features=features)
        ax, _ = record.plot(figure_width=15, figure_height=15)
        patch0 = mpatches.Patch(color='#ffd700', label='Good alleles')
        patch1 = mpatches.Patch(color='#ccccff', label='Other alleles')
        legend = ax.legend(handles=[patch0, patch1], loc='upper right', fontsize=10,
                           title=str(inputfile).replace('_', ' '),
                           title_fontsize=15)
        legend.get_frame().set_alpha(None)
        print(
            'Please introduce the name of the file for your image and extension')

        outputfile = input()
        ax.figure.savefig(str(outputfile), bbox_inches='tight')
