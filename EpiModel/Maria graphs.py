import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

phi_values = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Miniproject/Aphids/Copy of phi values['
                           '17].xlsx', sheet_name='Sheet4', index_col=0)
phi_values = phi_values.transpose()

fig, ax = plt.subplots(1, 1)

sns.lineplot(phi_values)
plt.axvline(x=0.008, color='#FF6EB4', label='Plateau threshold')
plt.xlabel('PhiA values', fontweight='bold', fontsize=12)
plt.ylabel('Percentage of cabbage infection (%)', fontweight='bold', fontsize=12)

plt.savefig('Optimization of Phi values.pdf', dpi=300)
plt.show()

# %%

aphids = pd.read_excel('/Users/u2176312/OneDrive - University of Warwick/Miniproject/Aphids/Copy of phi values['
                           '17].xlsx', sheet_name='Sheet2', index_col=0)


fig, ax = plt.subplots(1, 1)

sns.lineplot(aphids)
plt.axvline(x=200, color='#FF6EB4', label='Plateau threshold')
plt.xlabel('Number of viruliferous aphids', fontweight='bold', fontsize=12)
plt.ylabel('Percentage of cabbage infection (%)', fontweight='bold', fontsize=12)
plt.savefig('Optimization of viruliferous aphids.pdf', dpi=300)
plt.show()