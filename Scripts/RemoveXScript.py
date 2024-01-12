from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import argparse
import textwrap

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent('''\
    Creation of all peptides for CSP based on NCBI database.
    ------------------------------------------
    '''))

parser.add_argument("-i", "--in_file", metavar='file.fasta', dest="infile", help="Sequences in fasta", type=str)
parser.add_argument("-o", "--out_file", dest="outfile", help="Output file", type=str)


args = parser.parse_args()

fasta_file = SeqIO.parse(args.infile, 'fasta')

with open(args.outfile, 'a') as f:

    for seq_record in tqdm(fasta_file):
        sequence_id = seq_record.id
        sequence = seq_record.seq
        if ('X' in sequence) or ('Z' in sequence) or ('B' in sequence):
            continue
        else:
            new_seq_record = SeqRecord(sequence, id=sequence_id, description='')
            SeqIO.write(new_seq_record, f, 'fasta')
