# %% Combination of SNPs in the sequence and the MHC binding region
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def main():
    vcf_df = pd.read_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/SNP_INDEL_Pf3D7_03_v3.csv', index_col=0)
    P_falcifarum = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                               'CSP/Pfalciparum.genome.fasta', 'fasta')
    section_to_study = vcf_df[vcf_df['Start'] >= (81 * 3 + 221323)]  # Filters results to section of interest
    section_to_study = section_to_study[section_to_study['Stop'] <= (133 * 3 + 221323)]
    section_to_study = section_to_study[section_to_study['Effect'] != 'SYNONYMOUS_CODING']

    chromosome = []
    for seq_record in P_falcifarum:
        if seq_record.id == 'Pf3D7_03_v3':
            chromosome = seq_record.seq  # Get the sequence of the chromosome of interest
            break

    sequence_frame = chromosome[221323: 222516]  # Get the sequence of the region of interest

    relative_position = [x - 221323 for x in list(section_to_study['Pos'])]
    SNPs = list(section_to_study['Alt'].apply(eval))
    dict_SNPs = dict(zip(relative_position, SNPs))  # Create a dictionary linking relative position and SNPs
    indexed_sequence = pd.Series(sequence_frame)

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
                for i in range(len(combinations)):
                    # Replace the nucleotide at the SNP position with the SNP nucleotide
                    new_seq = list(combinations[i])
                    new_seq[pos] = nt
                    # Add the new combination to the list of new combinations
                    new_combinations.append("".join(new_seq))
            # Update the list of combinations with the new combinations
            combinations = new_combinations

        return combinations

    collector_object = []

    kmer = indexed_sequence  # Get the indexed kmer

    print('Generating all possible combinations...')
    new_sequence = generate_combinations(kmer.values, dict_SNPs)  # Generate the combination of sequences
    collector_object.extend(new_sequence)  # Add the new sequences to the list

    print('Creating fasta file of all possible combinations...')
    counter = 0
    with open('/Users/u2176312/OneDrive - University of Warwick/'
              'CSP/SNP_sequence_variation_nt.fasta', 'a') as f1:
        with open('/Users/u2176312/OneDrive - University of Warwick/'
                  'CSP/SNP_sequence_variation_aa.fasta', 'a') as f2:
            for sequence in tqdm(collector_object):
                nt_sequence = Seq(sequence)
                rv_complement = nt_sequence.reverse_complement()
                aa_sequence = rv_complement.translate(table=1)

                seq_record = SeqRecord(nt_sequence, id='Pf3D7_csp' + str(counter), description='')
                SeqIO.write(seq_record, f1, 'fasta')

                seq_record = SeqRecord(aa_sequence, id='Pf3D7_csp' + str(counter), description='')
                SeqIO.write(seq_record, f2, 'fasta')
                counter += 1


if __name__ == '__main__':
    main()
