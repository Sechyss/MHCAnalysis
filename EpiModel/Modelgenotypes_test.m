% Testing the new function of the simplified model

Population = 10000;
Ratio12 = 0.15;
Ratio11 = 0.15;
Ratio22 = 0.3;
Ratio33 = 0.4;
ratios = [Ratio12 Ratio11 Ratio22 Ratio33];

betavalues = [0.000015 0.000015 0.000015 0.000015];
time_new = 300;  % Time of introduction of the new strain 
Time_simulation = 365*10;  % Total time for the simulation (in days)
new_Strain = '1a2b';
seeds = 12;
parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

[out, t] = MultiGenotype_strain_model_2(Population, ratios, parameters);

saveas(gcf, '/Users/u2176312/OneDrive - University of Warwick/Model/New Model/Betatesting/Combination_newModel_noexternalinfection/Testing_Genotype_ratios.pdf');

