import pandas as pd
import pickle
import argparse
import textwrap
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
This script takes a multiple protein alignment file and a list of column names,
and returns a fasta file with kmer sequences with alignment coordinates based on a given reference ID.

Additionally it returns a dictionary with a key for each kmer and the IDs associated with that kmer.

'''))

parser.add_argument("-a", "--alignment", metavar='MUSCLE', dest="sequences", help="Alignment in fasta",
                    type=str, required=True)
parser.add_argument('-s', '--starting', metavar='starting_aa', dest="starting_point",
                    help="Starting aa position", type=int)
parser.add_argument("-r", "--reference", dest="ref_id", help="Reference ID", type=str,  required=True)
parser.add_argument("-o", "--output", metavar='file.fasta', dest="output", help="Output file to save",
                    type=str,  required=True)


args = parser.parse_args()

# #######################################################################################################################

alignment = AlignIO.read(args.sequences, format="fasta")
counter = args.starting_point
coordinates = []
for record in alignment:
    if record.id == args.ref_id:
        sequence = record.seq
        for aminoacid in sequence:
            if aminoacid != '-':
                coordinates.append(counter)
                counter += 1
            else:
                coordinates.append(counter)
    else:
        continue

df = pd.DataFrame(index=[records.id for records in alignment], columns=coordinates)

print('Collecting coordinates for ' + args.output)
for i, col in tqdm(enumerate(alignment)):
    aligment_id = col.id
    all_seq = col.seq
    for index, aminoacid in enumerate(all_seq):
        df.loc[aligment_id, coordinates[index]] = aminoacid

Gene_depth = {}
print('Estimating the gene depth for ' + args.output)
ending_pos = 0
for position in tqdm(range(len(df.columns))):
    if ending_pos != df.columns[-1]:
        starting_pos = df.columns[position]
        ending_pos = df.columns[position + 10]
        filtered_df = df[df.columns[position: position + 10]]
        count_rows_with_hyphen = filtered_df.eq('-').any(axis=1).sum()
        depth = len(df.index) - count_rows_with_hyphen
        Gene_depth.update({starting_pos: depth})

Sequences_dictionary = {}

print('Writing Kmer sequences for ' + args.output)
with open(str(args.output) + '_kmers.fasta', 'a') as f:

    for index, row in tqdm(df.iterrows()):
        filtered_row = row[row != '-']
        all_coordinates = filtered_row.index.tolist()
        full_sequence = "".join(filtered_row.values.tolist())

        ending_pos = 0
        for position in range(len(all_coordinates)):
            if ending_pos != all_coordinates[-1]:
                starting_pos = all_coordinates[position]
                ending_pos = all_coordinates[position + 10]

                kmer = "".join(filtered_row.values[position: position + 11])
                id_kmer = "".join([args.output, '_', str(starting_pos), '_', str(ending_pos), '_', index])
                if kmer not in Sequences_dictionary.keys():
                    Sequences_dictionary.update({kmer: [id_kmer]})
                else:
                    Sequences_dictionary[kmer].append(id_kmer)
                seq_record = SeqRecord(seq=Seq(kmer), id=id_kmer, description='')
                SeqIO.write(seq_record, f, 'fasta')
            else:
                continue


with open(str(args.output) + '.pickle', 'wb') as f2:
    pickle.dump(Sequences_dictionary, f2)

with open(str(args.output) + '_Kmer_depth.pickle', 'wb') as f3:
    pickle.dump(Gene_depth, f3)

with open(str(args.output) + '_kmers_filtered.fasta', 'a') as f4:
    Sequence_counter = 0

    for key in tqdm(Sequences_dictionary.keys()):
        peptide = key
        sequence_id = 'Sequence_' + str(Sequence_counter)
        seq_record = SeqRecord(Seq(peptide), id=sequence_id, description='')
        SeqIO.write(seq_record, f4, 'fasta')
        Sequence_counter += 1


