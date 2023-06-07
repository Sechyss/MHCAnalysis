import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%

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

# %% Plot log equation

from scipy import optimize


x = [100, 19015, 983, 664, 2, 58, 9179, 26588, 295, 31587, 503, 16389]
y = [91.9, 97.2, 66.6, 95.6, 8.1, 19.8, 40.2, 98.8, 26.5, 91.5, 98.0, 99.7]
x.sort()
y.sort()

def log_equation(x_value, intercept, slope):
    import numpy as np
    return intercept + slope * np.log(x_value)


p0 = [0.0001, 0.0001]
params, cv = optimize.curve_fit(log_equation, xdata=x, ydata=y, p0=p0)
plt.plot(x, y, 'o')
plt.plot(x, log_equation(x, params[0], params[1]), 'red')
plt.text(5000, 10, u'y = '+str(round(params[0],2))+'+'+str(round(params[1], 2))+'log(x)', style='italic')

plt.show()

