% Testing the new function of the simplified model

Population = 10000;
Ratio12 = 0;
Ratio11 = 0;
Ratio22 = 0;
Ratio33 = 1;
ratios = [Ratio12 Ratio11 Ratio22 Ratio33];

betavalues = [0.000015 0.000015 0.000015 0.000015];
time_new = 300;  % Time of introduction of the new strain 
Time_simulation = 365*3;  % Total time for the simulation (in days)
new_Strain = '1b2b';
seeds = 12;
parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

[out, t] = MultiGenotype_strain_model(Population, ratios, parameters);

