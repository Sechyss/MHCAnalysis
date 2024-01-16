%% Preparation of parameters

seeds = 1;

rng(seeds);

% Select the two genotypes to compare
Genotype1 = 2;
Genotype2 = 3;

proportions = zeros (1, 4);
for i = 1:10
    % Initialize randomNumbers and create an array with two random numbers
    randomNumbers = zeros(1, 4);
    firstrandomgenotype = rand(1);
    secondrandomgenotype = 1- firstrandomgenotype;
    randomNumbers(Genotype1) = firstrandomgenotype;
    randomNumbers(Genotype2) = secondrandomgenotype;
    
    % Create the vector
    proportions(i, :) = randomNumbers;
end

Population = 10000;
betavalues = [0.000015 0.000015 0.000015 0.000015];
time_new = 200;
Time_simulation = 365*3;
new_Strain = '1a2b';

parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

peaks = zeros(1, 4);

%% Running of all models with the different proportions

for i = 1:size(proportions, 1)
    clf;
    ratios = proportions(i, :);
    [out, t] = MultiGenotype_strain_model(Population, ratios, {betavalues, time_new, Time_simulation, new_Strain i});
    sum1a2a = max(out(2, :)+ out(14,:)+ out(20, :)+ out(28,:)+ out(34, :)+ out(42,:)+out(48, :));
    sum1a2b = max(out(3, :)+ out(15,:)+ out(21, :)+ out(29,:)+ out(35, :)+ out(43,:)+out(49, :));
    sum1b2a = max(out(4, :)+ out(16,:)+ out(22, :)+ out(30,:)+ out(36, :)+ out(44,:)+out(50, :));
    sum1b2b = max(out(5, :)+ out(17,:)+ out(23, :)+ out(31,:)+ out(37, :)+ out(45,:)+out(51, :));
    maxpeaks = [sum1a2a, sum1a2b, sum1b2a, sum1b2b];
    peaks = cat(1, peaks, maxpeaks);
end


peaks(1, :) = [];

%% Plotting of the increase of the population

figure;
hold on;

% Select the column of one of the genotypes to study 1 is Genotype 12, 
% 2 is 11, 3 is 22, and 4 is 33.
scatter(proportions(:,2), peaks(:,1), 50, 'red' , 'o', 'DisplayName','1a2a');
scatter(proportions(:,2), peaks(:,2), 50, 'green', 'o', 'DisplayName','1a2b');
scatter(proportions(:,2), peaks(:,3), 50, 'blue', 'o', 'DisplayName','1b2a');
scatter(proportions(:,2), peaks(:,4), 50, 'magenta', 'o', 'DisplayName','1b2b');


% Label the xlabel based on the genotypes, 1 is Genotype 12, 2 is 11, 3 is
% 22, and 4 is 33.
xlabel('Ratio of Genotype' )
ylabel('Max peak of infection')
legend('show')

hold off

