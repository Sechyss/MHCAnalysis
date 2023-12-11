% Testing the new function of the simplified model

Population = 10000;
Ratio12 = 0.1;
Ratio11 = 0.3;
Ratio22 = 0.3;
Ratio33 = 0.3;
ratios = [Ratio12 Ratio11 Ratio22 Ratio33];

betavalues = [0.002 0.002 0.002 0.002];
time_new = 200;
Time_simulation = 400;
new_Strain = '1b2b';
seeds = 42;
parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

[out, t] = MultiGenotype_strain_model(Population, ratios, parameters);

maxpeaks = max(out, [],2);
