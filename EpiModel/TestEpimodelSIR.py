import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, Nc.
Nc = 48
Na = 10
days_to_develop_infection_cabbage = 21

# Initial number of infected and recovered individuals, Ic0 and R0.
Ec0 = 0
Ic0 = 0
Ia0 = 10
# Everyone else, Sc0, is susceptible to infection initially.
Sc0 = Nc - Ic0
Sa0 = Na - Ia0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
phiC = 0.006
phiA = 0.004
alphaC = 1/days_to_develop_infection_cabbage

# A grid of time points (in days)
t = np.linspace(0, 44, 44)
# Rate of death and birth for the aphid population
b = 0.001
d = 0.001

# The SEI model differential equations.
def deriv(y, t, Ncabbage, Naphids, phiCabbage, phiAphids, birth, death, alpha):
    Sc, Sa, Ec, Ic, Ia = y
    dScdt = -(phiAphids * Ia) * Sc
    dEcdt = (phiAphids * Ia)* Sc - alpha*Ec
    dIcdt = alpha*Ec
    dSadt = -(death*phiCabbage*Ic) *Sa + birth*(Sa + Ia)
    dIadt = -death*Ia + (phiCabbage*Ic*Sa)
    return dScdt, dEcdt, dIcdt, dSadt, dIadt


# Initial conditions vector
y0 = Sc0, Sa0, Ec0, Ic0, Ia0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(Nc, Na, phiC, phiA, b, d, alphaC))
Scabbage, Ecabbage, Icabbage, Saphid, Iaphid = ret.T

# %% Plot the data on three separate curves for Scabagge(t), Icabbage(t), Saphid(t) and Iaphid(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, Scabbage, 'blue', alpha=0.5, lw=2, label='Susceptible Cabbages')
ax.plot(t, Ecabbage, 'purple', alpha=0.5, lw=2, label='Exposed Cabbages')
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
