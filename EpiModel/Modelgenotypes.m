% Testing the new function of the simplified model

Population = 10000;
Ratio12 = 0.4;
Ratio11 = 0.3;
Ratio22 = 0.2;
Ratio33 = 0.1;
ratios = [Ratio12 Ratio11 Ratio22 Ratio33];

betavalues = [0.002 0.002 0.002 0.002];
time_new = 200;  % Time of introduction of the new strain 
Time_simulation = 800;  % Total time for the simulation (in days)
new_Strain = '1b2a';
seeds = 42;
parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

MultiGenotype_strain_model(Population, ratios, parameters)
