from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

fasta_file = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                         'Otherproteins/Muscle_alignment/Thrombos_Nterminal_full_kmer_filtered.fasta', 'fasta')

with open('/Users/u2176312/OneDrive - University of Warwick/'
          'Otherproteins/Muscle_alignment/Thrombos_Nterminal_full_kmer_filtered_corrected.fasta', 'a') as f:

    for seq_record in fasta_file:
        sequence_id = seq_record.id
        sequence = seq_record.seq
        if ('X' in sequence) or ('Z' in sequence) or ('B' in sequence):
            continue
        else:
            new_seq_record = SeqRecord(sequence, id=sequence_id, description='')
            SeqIO.write(new_seq_record, f, 'fasta')
