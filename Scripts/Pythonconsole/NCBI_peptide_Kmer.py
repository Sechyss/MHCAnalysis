import argparse
import textwrap

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
    Creation of all peptides for CSP based on NCBI database.
    ------------------------------------------
    '''))

    parser.add_argument("-f", "--fasta", metavar='file.fasta', dest="sequences", help="Sequences in fasta", type=str)
    parser.add_argument("-k", "--kmer", dest="kmer", help="Kmer length", type=int)
    parser.add_argument("-o", "--output", metavar='file.fasta', dest="output", help="Output file to save", type=str)

    args = parser.parse_args()

    ####################################################################################################################
    fasta_file = SeqIO.parse(args.sequences, 'fasta')

    with open(args.output, 'a') as f1:

        for seq_record in fasta_file:
            starting_point = 0  # Starting point of the peptide sequence
            end_point = int(args.kmer)  # Ending point of the peptide sequence
            sequence_id = seq_record.id
            sequence = seq_record.seq
            while end_point <= len(sequence):
                peptide = Seq(sequence[starting_point:end_point])
                new_fasta = 'Kmer_' + str(starting_point) + '_' + str(end_point) + '_' + str(sequence_id)
                seq_record = SeqRecord(peptide, id=new_fasta, description='')
                SeqIO.write(seq_record, f1, 'fasta')
