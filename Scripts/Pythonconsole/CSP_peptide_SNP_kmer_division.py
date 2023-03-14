# /usr/bin/env python3
import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from MHCPipeline import generate_combinations

vcf_df = pd.read_csv('/Users/u2176312/OneDrive - University of Warwick/CSP/SNP_INDEL_Pf3D7_03_v3.csv', index_col=0)
P_falcifarum = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                           'CSP/Pfalciparum.genome.fasta', 'fasta')

locations = [int(re.findall(r'\d+', x)[0]) for x in vcf_df['Amino_acid_change']]

vcf_df['Amino_acid_change_location'] = locations

section_to_study = vcf_df[vcf_df['Effect'] != 'SYNONYMOUS_CODING']
section_to_study = section_to_study[section_to_study['Type'] == 'SNP']
section_to_study = section_to_study[
    section_to_study['Amino_acid_change_location'] >= 311]  # Filters results to region
section_to_study = section_to_study[section_to_study['Amino_acid_change_location'] <= 342]

chromosome = []
for seq_record in P_falcifarum:
    if seq_record.id == 'Pf3D7_03_v3':
        chromosome = seq_record.seq  # Get the sequence of the chromosome of interest
        break

sequence_frame = chromosome[221322: 222516]  # Get the sequence of the region of interest
reference_ordered = sequence_frame.reverse_complement()  # Get the reference sequence of the region of interest
reference_peptide = reference_ordered.translate(table=1)  # Get the reference sequence of the region of interest
reference_peptide = reference_peptide[311:342]

relative_aa_position = [int(x) for x in list(section_to_study['Amino_acid_change_location'])]

aminoacidchanged = [re.findall(r'\D+', x)[1] for x in section_to_study['Amino_acid_change']]
dict_SNPs = {}
for i in range(len(aminoacidchanged)):
    if relative_aa_position[i] not in dict_SNPs.keys():
        dict_SNPs.update({relative_aa_position[i]: [aminoacidchanged[i]]})
    else:
        dict_SNPs[relative_aa_position[i]].append(aminoacidchanged[i])

relative_aa_position = [x - 312 for x in dict_SNPs.keys()]
starting_point = 0  # Starting point of the MHC binding region
end_point = 10  # Ending point of the MHC binding region
indexed_sequence = pd.Series(reference_peptide)

collector_object = []
while end_point <= len(reference_peptide):  # Iterate through the region of interest
    kmer = indexed_sequence[starting_point: end_point]  # Get the indexed kmer
    SNP_hits = list(set(kmer.index) & set(relative_aa_position))  # Get the SNP hits in this kmer region
    SNP_hits.sort()  # Sort the SNP hits
    if len(SNP_hits) > 0:
        dict_sliced = {x - starting_point: dict_SNPs[x + 312] for x in
                       SNP_hits}  # Get the dictionary of matching regions
        new_sequence = generate_combinations(kmer.values, dict_sliced)  # Generate the combination of sequences
        collector_object.extend(new_sequence)  # Add the new sequences to the list
        starting_point += 1
        end_point += 1
    else:
        collector_object.append("".join(kmer))  # Add the reference sequence to the list
    starting_point += 1
    end_point += 1

counter = 0

with open('/Users/u2176312/OneDrive - University of Warwick/'
          'CSP/SNP_kmer_sequence_variation_311-342.fasta', 'a') as f1:
    for sequence in tqdm(collector_object):
        nt_sequence = Seq(sequence)

        if '*' not in nt_sequence:
            seq_record = SeqRecord(nt_sequence, id='Pf3D7_csp' + str(counter), description='')
            SeqIO.write(seq_record, f1, 'fasta')
            counter += 1
