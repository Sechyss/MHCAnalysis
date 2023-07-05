import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def delete_zero_sum_columns(df):
    """Deletes columns in a Pandas dataframe whose total sum is 0.

  Args:
    df: The Pandas dataframe.

  Returns:
    A new Pandas dataframe with the zero-sum columns deleted.
  """

    zero_sum_columns = []
    for column in df.columns:
        if df[column].sum() == 0:
            zero_sum_columns.append(column)

    new_df = df.drop(columns=zero_sum_columns)
    return new_df


diff = pd.read_excel('/Users/u2176312/OneDrive - University of '
                     'Warwick/CSP/Kmer_Variant_recognition_world_Comparison.xlsx', sheet_name='Diff')
diff = diff.drop(['Variants'], axis=1)

diff = delete_zero_sum_columns(diff)

diff = pd.DataFrame(diff).set_index('Kmer')
# Get the sum of each row
row_sums = diff.sum(axis=1)

# Delete the rows where the sum is 0
diff = diff[row_sums != 0]
diff = diff.reset_index()

df_melted = diff.melt(id_vars='Kmer', var_name='Countries', value_name='% of diff between top 3 & 2')

fig, ax = plt.subplots(figsize=(15, 10))

sns.barplot(x='Kmer', y='% of diff between top 3 & 2', hue='Countries', data=df_melted)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=9, fontweight='bold')
plt.yticks(np.arange(0, 110, 10))
plt.ylabel('% of diff between top 3 & 2', fontweight='bold', fontsize=15)
plt.yticks(fontsize=12, fontweight='bold')
plt.xlabel('List of Kmer peptides', fontweight='bold', fontsize=15)
ax.legend(loc='upper right', fontsize=15,
          title='Countries that changed',
          title_fontsize=20)
plt.tight_layout()
plt.savefig('/Users/u2176312/OneDrive - University of '
            'Warwick/CSP/Kmer_Variant_recognition_world_Comparison.pdf', dpi=300)
plt.show()
