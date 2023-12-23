% Testing the new function of the simplified model

Population = 10000;
Ratio12 = 0.5;
Ratio11 = 0.3;
Ratio22 = 0.2;
Ratio33 = 0.4;
ratios = [Ratio12 Ratio11 Ratio22 Ratio33];

betavalues = [0.00004 0.00008 0.00008 0.00008];
time_new = 200;  % Time of introduction of the new strain 
Time_simulation = 700;  % Total time for the simulation (in days)
new_Strain = '1b2b';
seeds = 123;
parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

[out, t] = MultiGenotype_strain_model(Population, ratios, parameters);
