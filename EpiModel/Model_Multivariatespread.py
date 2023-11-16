import numpy as np
import pandas as pd
from scipy.integrate import odeint


def sum_nested_tuple(nt):
    return sum(sum_nested_tuple(i) if isinstance(i, tuple) else i for i in nt)


# Selection of total Population genotype 12

S12_0 = 2600
I12_1a2a_0 = 15
I12_1a2b_0 = 15
I12_1b2a_0 = 15
I12_1b2b_0 = 15
J12_1a2a_0 = 0
J12_1a2b_0 = 0
J12_1b2a_0 = 0
J12_1b2b_0 = 0
R12_1a2a_0 = 0
R12_1a2b_0 = 0
R12_1b2a_0 = 0
R12_1b2b_0 = 0
H12_1a2a_0 = 0
H12_1a2b_0 = 0
H12_1b2a_0 = 0
H12_1b2b_0 = 0
M12_0 = 0

y_12 = (S12_0, I12_1b2a_0, I12_1a2a_0, I12_1a2b_0, I12_1b2b_0, J12_1a2a_0, J12_1a2b_0, J12_1b2a_0, J12_1b2b_0,
        R12_1a2a_0, R12_1a2b_0, R12_1b2a_0, R12_1b2b_0, H12_1a2a_0, H12_1a2b_0, H12_1b2a_0, H12_1b2b_0, M12_0)

# Selection of total Population genotype 11

S11_0 = 2600
I11_1a2a_0 = 15
I11_1a2b_0 = 15
I11_1b2a_0 = 15
I11_1b2b_0 = 15
J11_1a2a_0 = 0
J11_1a2b_0 = 0
J11_1b2a_0 = 0
J11_1b2b_0 = 0
R11_1a_0 = 0
R11_1b_0 = 0
H11_1a_0 = 0
H11_1b_0 = 0
M11_0 = 0

y_11 = (S11_0, I11_1b2a_0, I11_1a2a_0, I11_1a2b_0, I11_1b2b_0, J11_1a2a_0, J11_1a2b_0, J11_1b2a_0, J11_1b2b_0,
        R11_1a_0, R11_1b_0, H11_1a_0, H11_1b_0, M11_0)

# Selection of total Population genotype 22

S22_0 = 2600
I22_1a2a_0 = 15
I22_1a2b_0 = 15
I22_1b2a_0 = 15
I22_1b2b_0 = 15
J22_1a2a_0 = 0
J22_1a2b_0 = 0
J22_1b2a_0 = 0
J22_1b2b_0 = 0
R22_2a_0 = 0
R22_2b_0 = 0
H22_2a_0 = 0
H22_2b_0 = 0
M22_0 = 0

y_22 = (S22_0, I22_1b2a_0, I22_1a2a_0, I22_1a2b_0, I22_1b2b_0, J22_1a2a_0, J22_1a2b_0, J22_1b2a_0, J22_1b2b_0,
        R22_2a_0, R22_2b_0, H22_2a_0, H22_2b_0, M22_0)

# Selection of total Population genotype 22

S33_0 = 2600
I33_1a2a_0 = 15
I33_1a2b_0 = 15
I33_1b2a_0 = 15
I33_1b2b_0 = 15
M33_0 = 0

y_33 = (S33_0, I33_1b2a_0, I33_1a2a_0, I33_1a2b_0, I33_1b2b_0, M33_0)

# Creation of parameters

gamma_value = 1 / 7  # Gamma value represents force of recovery from infection in these case dependent on time (7 days)
sigma_value = 1 / 150  # Sigma value represents force of losing immunity against injection in these case dependent on
# time
birth_rate = 0.02  # Birthrate depending on adult population
death_rate = 0.002  # Death rate depending on external causes

beta_1a2a = 0.25  # Probability of infection by strain 1a2a
beta_1a2b = 0.15  # Probability of infection by strain 1a2b
beta_1b2a = 0.30  # Probability of infection by strain 1b2a
beta_1b2b = 0.11  # Probability of infection by strain 1b2b

beta_values = beta_1a2a, beta_1a2b, beta_1b2a, beta_1b2b

t = np.linspace(0, 800, 800)  # Time series for the simulation

y_values = y_12 + y_11 + y_22 + y_33  # This variable packs all the compartments into one

total_pop = sum_nested_tuple(y_values)


def deriv_equations(y: tuple, t, betas: tuple, gamma, sigma, death, birth, tap):
    # Unpackaging parameters from the input populations
    (S12, I12_1a2a, I12_1b2b, I12_1a2b, I12_1b2a, J12_1a2a, J12_1a2b, J12_1b2a, J12_1b2b, R12_1a2a, R12_1a2b, R12_1b2a,
     R12_1b2b, H12_1a2a, H12_1a2b, H12_1b2a, H12_1b2b, M12) = y[0:18]
    (S11, I11_1a2a, I11_1b2b, I11_1a2b, I11_1b2a, J11_1a2a, J11_1a2b, J11_1b2a, J11_1b2b, R11_1a, R11_1b, H11_1a, H11_1b
     , M11) = y[18:32]
    (S22, I22_1a2a, I22_1b2b, I22_1a2b, I22_1b2a, J22_1a2a, J22_1a2b, J22_1b2a, J22_1b2b, R22_2a, R22_2b, H22_2a, H22_2b
     , M22) = y[32:46]
    S33, I33_1a2a, I33_1b2b, I33_1a2b, I33_1b2a, M33 = y[46:]

    beta1a2a, beta1a2b, beta1b2a, beta1b2b = betas

    lambda1a2a = beta1a2a * (I12_1a2a + I11_1a2a + I22_1a2a + I33_1a2a + J12_1a2a + J11_1a2a + J22_1a2a)

    lambda1a2b = beta1a2b * (I12_1a2b + I11_1a2b + I22_1a2b + I33_1a2a + J12_1a2b + J11_1a2b + J22_1a2b)

    lambda1b2a = beta1b2a * (I12_1b2a + I11_1b2a + I22_1b2a + I33_1b2a + J12_1b2a + J11_1b2a + J22_1b2a)

    lambda1b2b = beta1b2b * (I12_1b2b + I11_1b2b + I22_1b2b + I33_1b2a + J12_1b2b + J11_1b2b + J22_1b2b)

    # Equations for genotype 12
    dS12dt = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S12 + birth * tap
    dI12_1a2a_dt = -(death + gamma) * I12_1a2a + lambda1a2a * S12
    dI12_1a2b_dt = -(death + gamma) * I12_1a2b + lambda1a2b * S12
    dI12_1b2a_dt = -(death + gamma) * I12_1b2a + lambda1b2a * S12
    dI12_1b2b_dt = -(death + gamma) * I12_1b2b + lambda1b2b * S12
    dR12_1a2a_dt = -(death + sigma) * R12_1a2a + gamma * I12_1a2a
    dR12_1a2b_dt = -(death + sigma) * R12_1a2b + gamma * I12_1a2b
    dR12_1b2a_dt = -(death + sigma) * R12_1b2a + gamma * I12_1b2a
    dR12_1b2b_dt = -(death + sigma) * R12_1b2b + gamma * I12_1b2b
    dH12_1a2a_dt = -(death + lambda1a2a) * H12_1a2a + sigma * R12_1a2a
    dH12_1a2b_dt = -(death + lambda1a2b) * H12_1a2a + sigma * R12_1a2b
    dH12_1b2a_dt = -(death + lambda1b2a) * H12_1b2a + sigma * R12_1b2a
    dH12_1b2b_dt = -(death + lambda1b2b) * H12_1b2b + sigma * R12_1b2b
    dJ12_1a2a_dt = -(death + gamma) * J12_1a2a + lambda1a2a * H12_1b2b
    dJ12_1a2b_dt = -(death + gamma) * J12_1a2b + lambda1a2b * H12_1b2a
    dJ12_1b2a_dt = -(death + gamma) * J12_1b2a + lambda1b2a * H12_1a2b
    dJ12_1b2b_dt = -(death + gamma) * J12_1b2b + lambda1b2b * H12_1a2a
    dM12_dt = -death * M12 + gamma * (I12_1a2a + I12_1b2a + I12_1a2b + I12_1b2b)

    # Equations for genotype 11 (or 13)
    dS11dt = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S11 + birth * tap
    dI11_1a2a_dt = -(death + gamma) * I11_1a2a + lambda1a2a * S11
    dI11_1a2b_dt = -(death + gamma) * I11_1a2b + lambda1a2b * S11
    dI11_1b2a_dt = -(death + gamma) * I11_1b2a + lambda1b2a * S11
    dI11_1b2b_dt = -(death + gamma) * I11_1b2b + lambda1b2b * S11
    dR11_1a_dt = -(death + sigma * R11_1a) + gamma * (I11_1a2a + I11_1a2b)
    dR11_1b_dt = -(death + sigma * R11_1b) + gamma * (I11_1b2a + I11_1b2b)
    dH11_1a_dt = -(death + lambda1a2a + lambda1a2b) * H11_1a + sigma * R11_1a
    dH11_1b_dt = -(death + lambda1b2a + lambda1b2b) * H11_1a + sigma * R11_1b
    dJ11_1a2a_dt = -(death + gamma) * J11_1a2a + lambda1a2a * H11_1b
    dJ11_1a2b_dt = -(death + gamma) * J11_1a2b + lambda1a2b * H11_1b
    dJ11_1b2a_dt = -(death + gamma) * J11_1b2a + lambda1b2a * H11_1a
    dJ11_1b2b_dt = -(death + gamma) * J11_1b2b + lambda1b2b * H11_1a
    dM11_dt = -death * M11 + gamma * (I11_1a2a + I11_1b2a + I11_1a2b + I11_1b2b)

    # Equations for genotype 22 (or 32)
    dS22dt = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S22 + birth * tap
    dI22_1a2a_dt = -(death + gamma) * I22_1a2a + lambda1a2a * S22
    dI22_1a2b_dt = -(death + gamma) * I22_1a2b + lambda1a2b * S22
    dI22_1b2a_dt = -(death + gamma) * I22_1b2a + lambda1b2a * S22
    dI22_1b2b_dt = -(death + gamma) * I22_1b2b + lambda1b2b * S22
    dR22_2a_dt = -(death + sigma * R22_2a) + gamma * (I22_1a2a + I22_1b2a)
    dR22_2b_dt = -(death + sigma * R22_2b) + gamma * (I22_1a2b + I22_1b2b)
    dH22_2a_dt = -(death + lambda1a2a + lambda1b2a) * H22_2a + sigma * R22_2a
    dH22_2b_dt = -(death + lambda1a2b + lambda1b2b) * H22_2a + sigma * R22_2b
    dJ22_1a2a_dt = -(death + gamma) * J22_1a2a + lambda1a2a * H22_2b
    dJ22_1a2b_dt = -(death + gamma) * J22_1a2b + lambda1a2b * H22_2a
    dJ22_1b2a_dt = -(death + gamma) * J22_1b2a + lambda1b2a * H22_2b
    dJ22_1b2b_dt = -(death + gamma) * J22_1b2b + lambda1b2b * H22_2a
    dM22_dt = -death * M22 + gamma * (I22_1a2a + I22_1b2a + I22_1a2b + I22_1b2b)

    # Equations for genotype 33
    dS33dt = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S33 + birth * tap + sigma * M33
    dI33_1a2a_dt = -(death + gamma) * I33_1a2a + lambda1a2a * S33
    dI33_1a2b_dt = -(death + gamma) * I33_1a2b + lambda1a2b * S33
    dI33_1b2a_dt = -(death + gamma) * I33_1b2a + lambda1b2a * S33
    dI33_1b2b_dt = -(death + gamma) * I33_1b2b + lambda1b2b * S33
    dM33_dt = -(death + sigma) * M33 + gamma * (I33_1a2a + I33_1b2a + I33_1a2b + I33_1b2b)

    return (dS12dt, dI12_1a2a_dt, dI12_1a2b_dt, dI12_1b2a_dt, dI12_1b2b_dt, dR12_1a2a_dt,
            dR12_1a2b_dt, dR12_1b2a_dt, dR12_1b2b_dt, dH12_1a2a_dt, dH12_1a2b_dt, dH12_1b2a_dt, dH12_1b2b_dt,
            dJ12_1a2a_dt, dJ12_1a2b_dt, dJ12_1b2a_dt, dJ12_1b2b_dt, dM12_dt, dS11dt, dI11_1a2a_dt, dI11_1a2b_dt,
            dI11_1b2a_dt, dI11_1b2b_dt, dR11_1a_dt, dR11_1b_dt, dH11_1a_dt, dH11_1b_dt, dJ11_1a2a_dt, dJ11_1a2b_dt,
            dJ11_1b2a_dt, dJ11_1b2b_dt, dM11_dt, dS22dt, dI22_1a2a_dt, dI22_1a2b_dt, dI22_1b2a_dt, dI22_1b2b_dt,
            dR22_2a_dt, dR22_2b_dt, dH22_2a_dt, dH22_2b_dt, dJ22_1a2a_dt, dJ22_1a2b_dt, dJ22_1b2a_dt, dJ22_1b2b_dt,
            dM22_dt, dS33dt, dI33_1a2a_dt, dI33_1a2b_dt, dI33_1b2a_dt, dI33_1b2b_dt, dM33_dt)


# ----------------------------------------------------------------
# In process to finish
ret1 = odeint(deriv_equations, y_values, t, args=(beta_values, gamma_value, sigma_value, death_rate, birth_rate,
                                                  total_pop))

dataframe = pd.DataFrame(data=ret1, columns=['dS12dt', 'dI12_1a2a_dt', 'dI12_1a2b_dt', 'dI12_1b2a_dt', 'dI12_1b2b_dt',
                                             'dR12_1a2a_dt', 'dR12_1a2b_dt', 'dR12_1b2a_dt', 'dR12_1b2b_dt',
                                             'dH12_1a2a_dt', 'dH12_1a2b_dt', 'dH12_1b2a_dt', 'dH12_1b2b_dt',
                                             'dJ12_1a2a_dt', 'dJ12_1a2b_dt', 'dJ12_1b2a_dt', 'dJ12_1b2b_dt', 'dM12_dt',
                                             'dS11dt', 'dI11_1a2a_dt', 'dI11_1a2b_dt', 'dI11_1b2a_dt', 'dI11_1b2b_dt',
                                             'dR11_1a_dt', 'dR11_1b_dt', 'dH11_1a_dt', 'dH11_1b_dt', 'dJ11_1a2a_dt',
                                             'dJ11_1a2b_dt', 'dJ11_1b2a_dt', 'dJ11_1b2b_dt', 'dM11_dt', 'dS22dt',
                                             'dI22_1a2a_dt', 'dI22_1a2b_dt', 'dI22_1b2a_dt', 'dI22_1b2b_dt',
                                             'dR22_2a_dt', 'dR22_2b_dt', 'dH22_2a_dt', 'dH22_2b_dt', 'dJ22_1a2a_dt',
                                             'dJ22_1a2b_dt', 'dJ22_1b2a_dt', 'dJ22_1b2b_dt', 'dM22_dt', 'dS33dt',
                                             'dI33_1a2a_dt', 'dI33_1a2b_dt', 'dI33_1b2a_dt', 'dI33_1b2b_dt', 'dM33_dt'])
