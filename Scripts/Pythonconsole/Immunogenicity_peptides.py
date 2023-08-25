from Bio import SeqIO
import pickle


def file_from_list(inputlist, outfile):
    with open(outfile, 'a') as file:
        for i in inputlist:
            file.write(str(i) + '\n')


fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                        'CSP/AllelePops/Kmer_CSP_region_273-375_aa.fasta', 'fasta')

listtofile = []
dict_tosave = {}
for seq_record in fastafile:
    listtofile.append(seq_record.seq)
    dict_tosave.update({seq_record.id: seq_record.seq})

file_from_list(listtofile, '/Users/u2176312/OneDrive - University of Warwick/'
                           'CSP/AllelePops/Immunogenicity_sequences_273-375.txt')
with open('/Users/u2176312/OneDrive - University of Warwick/'
          'CSP/AllelePops/Dictionary_ids_peptide_273-375.pickle', 'wb') as f:
    pickle.dump(dict_tosave, f)

# %% Repetition for the benign sequence from HLA ligand atlas

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                        'CSP/Humanrecogn/Immunnogenicity/hla_ligand_atlas/sequences_benign_human.fasta', 'fasta')

listtofile = []
dict_tosave = {}

for seq_record in fastafile:
    listtofile.append(seq_record.seq)
    dict_tosave.update({seq_record.id: seq_record.seq})

file_from_list(listtofile, '/Users/u2176312/OneDrive - University of Warwick/'
                           'CSP/Humanrecogn/Immunnogenicity/hla_ligand_atlas/sequences_benign_human.txt')
with open('/Users/u2176312/OneDrive - University of Warwick/'
          'CSP/Humanrecogn/Immunnogenicity/hla_ligand_atlas/sequences_benign_human.pickle', 'wb') as f:
    pickle.dump(dict_tosave, f)

# %% Repetition for the NCBI sequences

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/'
                        'CSP/NCBI_CSP/NCBI_CSP_peptides_11kmer_filtered.fasta', 'fasta')

listtofile = []
dict_tosave = {}

for seq_record in fastafile:
    listtofile.append(seq_record.seq)
    dict_tosave.update({seq_record.id: seq_record.seq})

file_from_list(listtofile, '/Users/u2176312/OneDrive - University of Warwick/'
                           'CSP/NCBI_CSP/Immunogenicity/NCBI_sequences_filtered.txt')
with open('/Users/u2176312/OneDrive - University of Warwick/'
          'CSP/NCBI_CSP/Immunogenicity/NCBI_sequences_filtered.pickle', 'wb') as f:
    pickle.dump(dict_tosave, f)
    