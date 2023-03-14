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
# Parse command line arguments
# -i FILE -t TASK -h HELP -o OUTPUT -th Threshold

parser.add_argument("-v", "--vcf", metavar='file.vcf.gz', dest="vcf", help="VCF file", type=str)
parser.add_argument("-m", "--mhc", metavar='file.txt', dest="mhc", help="MHC file", type=str)
parser.add_argument("-s", "--sequence", dest="fasta", help="Fasta file", type=str)
parser.add_argument("-o", "--output", metavar='file.pdf', dest="output", help="Output file to save", type=str)

args = parser.parse_args()


# =============================================================================

def main():
    # Features includes the name of the sequence, the length of the sequence, the number of elements...
    features = [GraphicFeature(start=1, end=18 * 3, color="olivedrab", label='Signal'),
                GraphicFeature(start=93 * 3, end=97 * 3, color="brown", label='Region I'),
                GraphicFeature(start=105 * 3, end=272 * 3, color="skyblue", label='Tandem repeats'),
                GraphicFeature(start=326 * 3, end=342 * 3, color="crimson", label='Th2R region'),
                GraphicFeature(start=366 * 3, end=380 * 3, color='green', label='Th3R region'),
                GraphicFeature(start=375 * 3, end=397 * 3, color="purple", label='GPI-anchor')]

    # Open input, add FILTER header, and open output file
    vcf_reader = pysam.VariantFile(str(args.vcf), 'r')

    data = []
    # Fetch region .  Note that the coordinates
    # in the API call are zero-based and describe half-open intervals.
    for record in vcf_reader.fetch('Pf3D7_03_v3', 221323, 222516):
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

    df = pd.DataFrame(data=data, columns=['Chrom', 'Pos', 'Start', 'Stop', 'ID', 'Ref', 'Alt', 'Type',
                                          'Qual', 'GQ',
                                          'Allele_frequency', 'Read_Depth', 'Amino_acid_change', 'Codon_change',
                                          'Effect',
                                          'Impact'])

    # Addition of the MHC table to compare SNPs and MHC binding sites

    df_MHC = pd.read_table(str(args.mhc), sep='\t')
    df_MHC_ranked = df_MHC[df_MHC['rank'] <= 1]
    df_MHC_ranked = df_MHC_ranked.sort_values(by=['start'])

    # Run this is you want the regions of binding with HLAs
    data_u = [[0, 0]]
    for index, row in df_MHC_ranked.iterrows():
        coordinates = [df_MHC.loc[index]['start'], df_MHC.loc[index]['end']]
        test_coordinates = range(df_MHC.loc[index]['start'], df_MHC.loc[index]['end'])

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
    # Use the data to create the figure of interest
    csp_nucleotide = SeqIO.parse(str(args.fasta), 'fasta')

    for sequence in csp_nucleotide:

        for element in data_u:
            corrected_start = int(element[0])
            corrected_end = int(element[1])

            snps = 0

            for index2, row2 in df.iterrows():
                position = int(re.findall(r'\d+', df.loc[index2]['Amino_acid_change'])[0])
                if (corrected_start <= position <= corrected_end) and ('SNP' in df.loc[index2]['Type']):
                    snps += 1

            features.append(GraphicFeature(start=corrected_start, end=corrected_end, color='grey',
                                           label='Number of SNPs: ' + str(snps)))

        to_plot = GraphicRecord(sequence_length=len(sequence.seq), features=features)
        ax, _ = to_plot.plot()

    plt.savefig(str(args.output), dpi=300)
    plt.show()


if __name__ == '__main__':
    main()
