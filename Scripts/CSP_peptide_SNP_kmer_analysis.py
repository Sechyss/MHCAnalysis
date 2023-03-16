# /usr/bin/env python3
import argparse
import re
import textwrap

import pandas as pd
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
    Creation of all potential sequences for CSP based on SNPs.
    ------------------------------------------
    '''))
    # Parse command line arguments
    # -v vcftable -g Malariagenome -h HELP -o OUTPUT -r rangeMHC

    parser.add_argument("-v", "--vcf", metavar='file.csv', dest="vcftable", help="VCF table", type=str)
    parser.add_argument("-g", "--fasta", metavar='file.fasta', dest="genome", help="Genome in fasta", type=str)
    parser.add_argument("-s", "--start", dest="start", help="MHC region start", type=int)
    parser.add_argument("-e", "--end", dest="end", help="MHC region end", type=int)
    parser.add_argument("-o", "--output", metavar='file.pdf', dest="output", help="Output file to save", type=str)

    args = parser.parse_args()

    ####################################################################################################################
    def generate_combinations(ref_seq, snp_dict):
        """
        Generate all possible combinations of the reference sequence and SNPs dictionary.

        Parameters:
        ref_seq (pd.Series): Pandas Series representing the reference DNA sequence.
        snp_dict (dict): Dictionary with SNPs using the position as keys and a list of nucleotides as the value.

        Returns:
        list: List of all possible combinations as strings.
        """
        # Initialize the list of all possible combinations with the reference sequence
        combinations = ["".join(ref_seq)]

        # Loop through the SNP dictionary and generate all possible combinations
        for pos, nts in snp_dict.items():
            # Check that the SNP position is within the range of valid indices for the reference sequence
            if pos < 0 or pos >= len(ref_seq):
                raise ValueError(f"Invalid SNP position: {pos}")
            # Make a copy of the list of combinations
            new_combinations = combinations.copy()
            # Loop through the list of combinations and generate a new combination for each SNP
            for nt in nts:
                for y in range(len(combinations)):
                    # Replace the nucleotide at the SNP position with the SNP nucleotide
                    new_seq = list(combinations[y])
                    new_seq[pos] = nt
                    # Add the new combination to the list of new combinations
                    new_combinations.append("".join(new_seq))
            # Update the list of combinations with the new combinations
            combinations = new_combinations

        return combinations

    vcf_df = pd.read_csv(args.vcftable, index_col=0)
    P_falcifarum = SeqIO.parse(args.genome, 'fasta')
    mhc_region_start = int(args.start)
    mhc_region_end = int(args.end)

    locations = [int(re.findall(r'\d+', x)[0]) for x in vcf_df['Amino_acid_change']]

    vcf_df['Amino_acid_change_location'] = locations

    section_to_study = vcf_df[vcf_df['Effect'] != 'SYNONYMOUS_CODING']
    section_to_study = section_to_study[section_to_study['Type'] == 'SNP']
    section_to_study = section_to_study[
        section_to_study['Amino_acid_change_location'] >= mhc_region_start]  # Filters results to region
    section_to_study = section_to_study[section_to_study['Amino_acid_change_location'] <= mhc_region_end]

    chromosome = []
    for seq_record in P_falcifarum:
        if seq_record.id == 'Pf3D7_03_v3':
            chromosome = seq_record.seq  # Get the sequence of the chromosome of interest
            break

    sequence_frame = chromosome[221322: 222516]  # Get the sequence of the region of interest
    reference_ordered = sequence_frame.reverse_complement()  # Get the reference sequence of the region of interest
    reference_peptide = reference_ordered.translate(table=1)  # Get the reference sequence of the region of interest
    reference_peptide = reference_peptide[mhc_region_start:mhc_region_end]

    relative_aa_position = [int(x) for x in list(section_to_study['Amino_acid_change_location'])]

    aminoacidchanged = [re.findall(r'\D+', x)[1] for x in section_to_study['Amino_acid_change']]
    dict_SNPs = {}
    for i in range(len(aminoacidchanged)):
        if relative_aa_position[i] not in dict_SNPs.keys():
            dict_SNPs.update({relative_aa_position[i]: [aminoacidchanged[i]]})
        else:
            dict_SNPs[relative_aa_position[i]].append(aminoacidchanged[i])

    relative_aa_position = [x - int(mhc_region_start + 1) for x in dict_SNPs.keys()]
    starting_point = 0  # Starting point of the MHC binding region
    end_point = 11  # Ending point of the MHC binding region
    indexed_sequence = pd.Series(reference_peptide)

    collector_object = []
    collector_dict = {}
    while end_point <= len(reference_peptide):  # Iterate through the region of interest
        kmer = indexed_sequence[starting_point: end_point]  # Get the indexed kmer
        SNP_hits = list(set(kmer.index) & set(relative_aa_position))  # Get the SNP hits in this kmer region
        SNP_hits.sort()  # Sort the SNP hits
        if len(SNP_hits) > 0:
            dict_sliced = {x - starting_point: dict_SNPs[x + int(mhc_region_start + 1)] for x in
                           SNP_hits}  # Get the dictionary of matching regions
            new_sequence = generate_combinations(kmer.values, dict_sliced)  # Generate the combination of sequences
            collector_object.extend(new_sequence)  # Add the new sequences to the list
            collector_dict.update(
                {'Kmer' + str(mhc_region_start + starting_point) + '-' + str(mhc_region_start + end_point):
                     new_sequence})
            starting_point += 1
            end_point += 1
        else:
            collector_object.append("".join(kmer))  # Add the reference sequence to the list
            collector_dict.update(
                {'Kmer' + str(mhc_region_start + starting_point) + '-' + str(mhc_region_start + end_point):
                     "".join(kmer)})
            starting_point += 1
            end_point += 1

    counter = 0

    with open(args.output, 'a') as f1:  # Open file to write the Fasta sequences
        for sequence in tqdm(collector_object):
            nt_sequence = Seq(sequence)

            if '*' not in nt_sequence:  # Remove any incomplete sequence due to stop codons
                seq_record = SeqRecord(nt_sequence, id='Pf3D7_csp' + str(counter), description='')
                SeqIO.write(seq_record, f1, 'fasta')
                counter += 1
            else:
                counter += 1
    with open(args.output + '_dict.pickle', 'wb') as f2:
        pickle.dump(collector_dict, f2)  # Save the dictionary to a pickle file for future use


if __name__ == '__main__':
    main()
