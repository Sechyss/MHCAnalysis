import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, Nc.
Nc = 1000
Na = 1000
# Initial number of infected and recovered individuals, Ic0 and R0.
Ic0 = 0
Ia0 = 10
# Everyone else, Sc0, is susceptible to infection initially.
Sc0 = Nc - Ic0
Sa0 = Na - Ia0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
phiC = 0.02
phiA = 0.04
# A grid of time points (in days)
t = np.linspace(0, 365, 365)
# Rate of death and birth for the aphid population
b = 0.10
d = 0.08


# The SIR model differential equations.
def deriv(y, t, Ncabbage, Naphids, phiCabbage, phiAphids, birth, death):
    Sc, Sa, Ic, Ia = y
    dScdt = -(phiAphids * Sc * Ia) / Ncabbage
    dIcdt = (phiAphids * Sc * Ia) / Ncabbage
    dSadt = birth * (Sa+Ia) /Naphids - (phiCabbage * Ic + death) * Sa/Naphids
    dIadt = (phiCabbage * Ic) * Sa/Naphids - death*(Ia/Naphids)
    return dScdt, dIcdt, dSadt, dIadt


# Initial conditions vector
y0 = Sc0, Ic0, Sa0, Ia0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(Nc, Na, phiC, phiA, b, d))
Scabbage, Icabbage, Saphid, Iaphid = ret.T

# Plot the data on three separate curves for Scabagge(t), Icabbage(t), Saphid(t) and Iaphid(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, Scabbage, 'blue', alpha=0.5, lw=2, label='Susceptible Cabbages')
ax.plot(t, Icabbage, 'red', alpha=0.5, lw=2, label='Infected Cabbages')
ax.plot(t, Saphid, 'green', alpha=0.5, lw=2, label='Susceptible Aphids')
ax.plot(t, Iaphid, 'black', alpha=0.5, lw=2, label='Infected Aphids')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number')
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()
