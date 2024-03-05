#!/usr/bin/env python3

""" This python script produces a figure that summarizes the amount of SNPs in a particular protein
sequence into a single plot.

It is important to make sure that all the packages are downloaded before running this script and in the correct version,
otherwise errors may be raised by the computer.
"""

# Upload of packages to run the script
import argparse
import textwrap
import warnings

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import pysam
import re
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord

warnings.simplefilter(action='ignore', category=FutureWarning)
# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Analysis of MHC prediction results in combination with SNPs in a particular protein sequence.
------------------------------------------
'''))

parser.add_argument("-v", "--vcf", metavar='file.vcf.gz', dest="vcf", help="VCF file", type=str)
parser.add_argument("-m", "--mhc", metavar='file.txt', dest="mhc", help="MHC file", type=str)
parser.add_argument("-s", "--sequence", dest="fasta", help="Fasta file", type=str)
parser.add_argument("-o", "--output", metavar='file.pdf', dest="output", help="Output file to save", type=str)

args = parser.parse_args()


# =============================================================================

def main():
    # Features includes the name of the sequence, the length of the sequence, the number of elements...
    features = [GraphicFeature(start=1, end=18, color="olivedrab", label='Signal'),
                GraphicFeature(start=93, end=97, color="brown", label='Region I'),
                GraphicFeature(start=105, end=272, color="skyblue", label='Tandem repeats'),
                GraphicFeature(start=326, end=342, color="crimson", label='Th2R region'),
                GraphicFeature(start=366, end=380, color='green', label='Th3R region'),
                GraphicFeature(start=375, end=397, color="purple", label='GPI-anchor')]

    # Open input, add FILTER header, and open output file
    vcf_reader = pysam.VariantFile(str(args.vcf), 'r')

    data = []
    # Fetch region .  Note that the coordinates
    # in the API call are zero-based and describe half-open intervals.
    for record in vcf_reader.fetch('Pf3D7_03_v3', 221323, 222516):
        # Extract the information from the record and add it to the coordinates of features for this region
        if 'SNPEFF_GENE_NAME' in record.info.keys() and record.info['SNPEFF_GENE_NAME'] == 'CSP':
            if record.rlen == 1:
                SNP_type = 'SNP'
            else:
                SNP_type = 'INDEL/OVERLAP'
            if 'SNPEFF_CODON_CHANGE' in record.info.keys():
                codon_change = record.info['SNPEFF_CODON_CHANGE']
            else:
                codon_change = '-'
            data.append(
                [record.chrom, record.pos, record.start, record.stop, record.id, record.ref, list(record.alts),
                 SNP_type,
                 record.qual,
                 record.info['GC'], record.info['AF'], record.info['DP'],
                 record.info['SNPEFF_AMINO_ACID_CHANGE'],
                 codon_change, record.info['SNPEFF_EFFECT'],
                 record.info['SNPEFF_IMPACT']])
    # Group the regions into a single DataFrame to create a CSV file if necessary, this can be added later
    df = pd.DataFrame(data=data, columns=['Chrom', 'Pos', 'Start', 'Stop', 'ID', 'Ref', 'Alt', 'Type',
                                          'Qual', 'GQ',
                                          'Allele_frequency', 'Read_Depth', 'Amino_acid_change', 'Codon_change',
                                          'Effect',
                                          'Impact'])

    # Addition of the MHC table from predictionbinding.py to compare SNPs and MHC binding sites

    df_MHC = pd.read_table(str(args.mhc), sep='\t')
    df_MHC_ranked = df_MHC[df_MHC['rank'] <= 1]
    df_MHC_ranked = df_MHC_ranked.sort_values(by=['start'])

    # Grouping of the MHC binding regions regardless of HLA that binds
    data_u = [[0, 0]]
    for index, row in df_MHC_ranked.iterrows():
        # Extract the coordinates of the SNPs and the potential points in between
        coordinates = [df_MHC.loc[index]['start'], df_MHC.loc[index]['end']]
        test_coordinates = range(df_MHC.loc[index]['start'], df_MHC.loc[index]['end'])

        # Counter is used to calculate if regions need to be further compared
        counter = 0
        temporary_addition_min = []
        temporary_addition_max = []
        for element in data_u:
            extracted_coordinates = set(range(int(element[0]), int(element[1])))

            if len(extracted_coordinates.intersection(test_coordinates)) != 0:
                counter += 1
                index = data_u.index(element)
                temporary_addition_min.append(min([int(element[0]), int(coordinates[1])]))
                temporary_addition_max.append(max([int(element[0]), int(coordinates[1])]))

            else:
                continue

        if counter == 0:
            data_u.append(coordinates)
        else:
            data_u[index] = [min(temporary_addition_min), max(temporary_addition_max)]
    data_u.pop(0)
    print('The regions of the protein that contain MHC binding regions are grouped as follows:')
    print(', '.join(map(str, data_u)))
    # Use the data to create the figure of interest
    csp_nucleotide = SeqIO.parse(str(args.fasta), 'fasta')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), height_ratios=[3, 1])

    features2 = []
    for sequence in csp_nucleotide:

        for element in data_u:
            corrected_start = int(element[0])
            corrected_end = int(element[1])
            snps = 0
            high_snps = 0
            moderate_snps = 0
            low_snps = 0
            for index2, row2 in df.iterrows():
                position = int(re.findall(r'\d+', df.loc[index2]['Amino_acid_change'])[0])
                if (corrected_start <= position <= corrected_end) and ('SNP' in df.loc[index2]['Type']):
                    snps += 1
                    if df.loc[index2]['Effect'] == 'NON_SYNONYMOUS_CODING':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='red'))
                        moderate_snps += 1
                    elif df.loc[index2]['Effect'] == 'CODON_INSERTION':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='navy'))
                        moderate_snps += 1
                    elif df.loc[index2]['Effect'] == 'CODON_DELETION':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='magenta'))
                    elif df.loc[index2]['Effect'] == 'CODON_CHANGE_PLUS_CODON_DELETION':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='yellow'))
                        moderate_snps += 1
                    elif df.loc[index2]['Effect'] == 'CODON_CHANGE_PLUS_CODON_INSERTION':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='orange'))
                        moderate_snps += 1
                    elif df.loc[index2]['Effect'] == 'FRAME_SHIFT':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='darkviolet'))
                        high_snps += 1
                    elif df.loc[index2]['Effect'] == 'STOP_GAINED':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='black'))
                        high_snps += 1
                    elif df.loc[index2]['Effect'] == 'SYNONYMOUS_CODING':
                        features2.append(
                            GraphicFeature(start=position,
                                           end=position + 1, color='forestgreen'))
                        low_snps += 1

            features.append(
                GraphicFeature(start=corrected_start,
                               end=corrected_end, color='grey'))

        to_plot = GraphicRecord(sequence_length=len(sequence.seq), features=list(set(features)))
        to_plot2 = GraphicRecord(sequence_length=len(sequence.seq), features=list(set(features2)))
        ax, _ = to_plot.plot(ax=ax1, with_ruler=False)
        ax2, _ = to_plot2.plot(ax=ax2, with_ruler=True)
        patch0 = mpatches.Patch(color='grey', label='MHC binding region')
        patch1 = mpatches.Patch(color='forestgreen', label='Synonymous coding')
        patch2 = mpatches.Patch(color='red', label='Non-synonymous coding')
        patch3 = mpatches.Patch(color='orange', label='Codon change plus codon insertion')
        patch4 = mpatches.Patch(color='yellow', label='Codon change plus codon deletion')
        patch5 = mpatches.Patch(color='navy', label='Codon insertion')
        patch6 = mpatches.Patch(color='magenta', label='Codon deletion')
        patch7 = mpatches.Patch(color='black', label='Stop gained')
        patch8 = mpatches.Patch(color='darkviolet', label='Frame shift')
        list_handles = [patch0, patch1, patch2, patch3, patch4, patch5, patch6, patch7, patch8]
        legend = fig.legend(handles=list_handles, loc='upper left')
        legend.get_frame().set_alpha(None)
        ax1.set_title('MHC binding region', loc='right', y=0.7, weight='bold')
        ax2.set_title('SNPs density', loc='right', y=0.55, weight='bold')
    print('Saving figure to file')
    plt.savefig(str(args.output), dpi=300)
    plt.show()


if __name__ == '__main__':
    main()
