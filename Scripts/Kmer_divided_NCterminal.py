import argparse
import textwrap
import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
    Creation of all peptides for proteins based on NCBI database.
    ------------------------------------------
    '''))

    parser.add_argument("-f", "--fastaN", metavar='fileN.fasta', dest="sequencesN", help="Sequences in N-terminal fasta",
                        type=str)
    parser.add_argument("-g", "--fastaC", metavar='fileC.fasta', dest="sequencesC", help="Sequences in C-terminal fasta",
                        type=str)
    parser.add_argument("-k", "--kmer", dest="kmer", help="Kmer length", type=int)
    parser.add_argument("-o", "--outputN", metavar='fileN.fasta', dest="outputN", help="Output file to save N-terminal",
                        type=str)
    parser.add_argument("-u", "--outputC", metavar='fileC.fasta', dest="outputC", help="Output file to save C-terminal",
                        type=str)
    args = parser.parse_args()

    fasta_file = SeqIO.parse(args.sequencesN, 'fasta')
    dictionary = {}
    with open(str(args.outputN) + '.fasta', 'a') as f1:

        for seq_record in tqdm(fasta_file):
            starting_point = 0  # Starting point of the peptide sequence
            end_point = int(args.kmer)  # Ending point of the peptide sequence
            sequence_id = seq_record.id
            sequence = seq_record.seq
            while end_point <= len(sequence):
                peptide = Seq(sequence[starting_point:end_point])
                new_fasta = 'Kmer_' + str(starting_point) + '_' + str(end_point - 1) + '_' + str(sequence_id)
                seq_record = SeqRecord(peptide, id=new_fasta, description='')
                SeqIO.write(seq_record, f1, 'fasta')
                if peptide in dictionary.keys():
                    dictionary[peptide].append(new_fasta)
                else:
                    dictionary.update({peptide: [new_fasta]})
                starting_point += 1
                end_point += 1

    with open(str(args.outputN) + '_filtered.fasta', 'a') as f2:
        Sequence_counter = 0

        for key in tqdm(dictionary.keys()):
            peptide = key
            sequence_id = 'Sequence_' + str(Sequence_counter)
            seq_record = SeqRecord(peptide, id=sequence_id, description='')
            SeqIO.write(seq_record, f2, 'fasta')
            Sequence_counter += 1

    with open(str(args.outputN) + '.pickle', 'wb') as f:
        pickle.dump(dictionary, f)

    fasta_file = SeqIO.parse(args.sequencesC, 'fasta')
    dictionary = {}

    with open(str(args.outputC) + '.fasta', 'a') as f1:

        for seq_record in tqdm(fasta_file):
            starting_point = 0  # Starting point of the peptide sequence
            end_point = int(args.kmer)  # Ending point of the peptide sequence
            sequence_id = seq_record.id
            sequence = seq_record.seq
            while end_point <= len(sequence):
                peptide = Seq(sequence[starting_point:end_point])
                new_fasta = 'Kmer_' + str(starting_point) + '_' + str(end_point - 1) + '_' + str(sequence_id)
                seq_record = SeqRecord(peptide, id=new_fasta, description='')
                SeqIO.write(seq_record, f1, 'fasta')
                if peptide in dictionary.keys():
                    dictionary[peptide].append(new_fasta)
                else:
                    dictionary.update({peptide: [new_fasta]})
                starting_point += 1
                end_point += 1

    with open(str(args.outputC) + '_filtered.fasta', 'a') as f2:
        for key in tqdm(dictionary.keys()):
            peptide = key
            sequence_id = 'Sequence_' + str(Sequence_counter)
            seq_record = SeqRecord(peptide, id=sequence_id, description='')
            SeqIO.write(seq_record, f2, 'fasta')
            Sequence_counter += 1

    with open(str(args.outputC) + '.pickle', 'wb') as f:
        pickle.dump(dictionary, f)


if __name__ == '__main__':
    main()
