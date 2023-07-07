from Bio import SeqIO

hypothetical_peptide_CSP = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                                       '/AllelePops/Kmer_CSP_region_273-375_aa.fasta', 'fasta')
real_peptide_CSP = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                               'NCBI_CSP/NCBI_CSP_peptides_11kmer_filtered.fasta', 'fasta')

raw_real_sequences = []
for seq_record in real_peptide_CSP:
    sequence = seq_record.seq
    raw_real_sequences.append(sequence)

raw_hypo_sequences = []
for seq_record in hypothetical_peptide_CSP:
    sequence = seq_record.seq
    raw_hypo_sequences.append(sequence)

matching_seqs = set(raw_real_sequences).intersection(set(raw_hypo_sequences))

print(str(round(len(matching_seqs)/len(raw_real_sequences), 2)) + '%  matching with the real sequences')

print(str(round(len(matching_seqs)/len(raw_hypo_sequences), 2)) + '%  matching with the hypothetical sequences')
