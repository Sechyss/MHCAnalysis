import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from dna_features_viewer import GraphicFeature, GraphicRecord
from Bio import SeqIO

Table = pd.read_table('/Users/u2176312/OneDrive - University of Warwick/CSP/'
                      'NCBI_CSP/NCBI_Pf3D7_blastp.tsv', sep='\t', header=None)

fastafile = SeqIO.parse('/Users/u2176312/OneDrive - University of Warwick/CSP/CSP_Pf3D7.faa', 'fasta')

features = [GraphicFeature(start=1, end=18, color="olivedrab", label='Signal'),
            GraphicFeature(start=93, end=97, color="brown", label='Region I'),
            GraphicFeature(start=105, end=272, color="skyblue", label='Tandem repeats'),
            GraphicFeature(start=322, end=375, color="crimson", label='TSP type-1'),
            GraphicFeature(start=375, end=397, color="purple", label='GPI-anchor')]
features2 = [GraphicFeature(start=1, end=18, color="olivedrab", label='Signal'),
             GraphicFeature(start=93, end=97, color="brown", label='Region I'),
             GraphicFeature(start=105, end=272, color="skyblue", label='Tandem repeats'),
             GraphicFeature(start=322, end=375, color="crimson", label='TSP type-1'),
             GraphicFeature(start=375, end=397, color="purple", label='GPI-anchor')]
sequence = []
for seq_record in fastafile:
    sequence = seq_record.seq

    starts = list(set(Table[6].tolist()))
    for item in starts:
        temptable = Table[Table[6] == item]
        endPoint = int(temptable[7].tolist()[0])
        counter = len(temptable[6].tolist())
        features.append(GraphicFeature(start=item, end=endPoint, color="#ffd700", label=str(counter)))

    Table2 = pd.read_csv('/Users/u2176312/OneDrive - University of Warwick/'
                         'CSP/Netchop/result_CSP_Pf3D7_netchop.csv', header=0)
    filteredTable2 = Table2[Table2['prediction_score'] >= 0.5]

    for index, row in filteredTable2.iterrows():
        startingposition = row['#']
        endposition = startingposition + 1
        features2.append(GraphicFeature(start=startingposition, end=endposition, color='#ccccff'))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), height_ratios=[3, 1])
to_plot = GraphicRecord(sequence_length=len(sequence), features=features)
to_plot2 = GraphicRecord(sequence_length=len(sequence), features=features2)
ax, _ = to_plot.plot(ax=ax1, with_ruler=False)
ax2, _ = to_plot2.plot(ax=ax2, with_ruler=True)
ax1.set_title('NCBI Blastp hits', loc='left', y=0.7, weight='bold')
ax2.set_title('NetChop hits', loc='left', y=0.55, weight='bold')

plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of Warwick/'
            'CSP/Netchop/Plot_Blastp_NCBI_Pf3D7_Netchop.pdf', dpi=600)
plt.show()

