%% Establish the parameters of the model
Population = 10000;
proportions = [0.25 0.25 0.25 0.25];
time_new = 200;
Time_simulation = 365*3;
new_Strain = '1b2a';

%% Creation of the different beta value combinations

seeds = 123;

rng(seeds);

% Initialize variables
numSamples = 1000;
betavalues = zeros(numSamples, 4);

% Define valid ranges for each element
validRanges = [0.0000075, 0.00007125];

for i = 1:numSamples
    % Generate random numbers within the valid ranges
    randomNumbers = rand(1, 4) .* diff(validRanges, 1, 2) + validRanges(:, 1)';
    
    % Store the valid combination
    betavalues(i, :) = randomNumbers;
end

%% Ratios of beta values
ratio1a2b = betavalues(:, 2) ./ betavalues (:, 1);
min_value = min(ratio1a2b);
max_value = max(ratio1a2b);

normalized_ratio1a2b = (ratio1a2b - min_value) / (max_value - min_value);


ratio1b2a = betavalues(:, 3) ./ betavalues (:, 1);
min_value = min(ratio1b2a);
max_value = max(ratio1b2a);

normalized_ratio1b2a = (ratio1b2a - min_value) / (max_value - min_value);

ratio1b2b = betavalues(:, 4) ./ betavalues (:, 1);
min_value = min(ratio1b2b);
max_value = max(ratio1b2b);

normalized_ratio1b2b = (ratio1b2b - min_value) / (max_value - min_value);


%% Run the model

peaks = zeros(1,  4);
for i = 1:size(betavalues, 1)
    clf;
    testingbetavalue = betavalues(i, :);
    [out, t] = MultiGenotype_strain_model(Population, proportions, {testingbetavalue, time_new, Time_simulation, new_Strain, seeds});
    sum1a2a = max(out(2, :)+ out(14,:)+ out(20, :)+ out(28,:)+ out(34, :)+ out(42,:)+out(48, :));
    sum1a2b = max(out(3, :)+ out(15,:)+ out(21, :)+ out(29,:)+ out(35, :)+ out(43,:)+out(49, :));
    sum1b2a = max(out(4, :)+ out(16,:)+ out(22, :)+ out(30,:)+ out(36, :)+ out(44,:)+out(50, :));
    sum1b2b = max(out(5, :)+ out(17,:)+ out(23, :)+ out(31,:)+ out(37, :)+ out(45,:)+out(51, :));
    maxpeaks = [sum1a2a, sum1a2b, sum1b2a, sum1b2b];
    peaks = cat(1, peaks, maxpeaks);
end

peaks(1, :) = [];
%% Plotting of the results from the model

figure;

hold on;
scatter(betavalues(:, 1), peaks(:,1), 50, 'red' , 'o', 'DisplayName','1a2a');
scatter(betavalues(:, 2), peaks(:,2), 50, 'green', 'o', 'DisplayName','1a2b');
scatter(betavalues(:, 3), peaks(:,3), 50, 'blue', 'o', 'DisplayName','1b2a');
scatter(betavalues(:, 4), peaks(:,4), 50, 'magenta', 'o', 'DisplayName','1b2b');

xlabel('Beta values');
ylabel('Max peak of infection');
legend('show');

hold off;

saveas(gcf, '/Users/u2176312/OneDrive - University of Warwick/Model/New Model/Betatesting/1b2a_infection_200days_combined_betavalues.pdf')

%% Plotting of the results of the peaks by the ratios of the beta values
    
figure;

hold on;
if strcmpi(new_Strain, '1a2b')
    scatter(normalized_ratio1a2b, peaks(:,1), 50, 'red' , 'o', 'DisplayName','1a2a');
    scatter(normalized_ratio1a2b, peaks(:,2), 50, 'green', 'o', 'DisplayName','1a2b');
elseif strcmp (new_Strain, '1b2a')
    scatter(normalized_ratio1b2a, peaks(:,1), 50, 'red' , 'o', 'DisplayName','1a2a');
    scatter(normalized_ratio1b2a, peaks(:,3), 50, 'blue', 'o', 'DisplayName','1b2a');
elseif strcmp (new_Strain, '1b2b')
    scatter(normalized_ratio1b2b, peaks(:,1), 50, 'red' , 'o', 'DisplayName','1a2a');
    scatter(normalized_ratio1b2b, peaks(:,4), 50, 'magenta', 'o', 'DisplayName','1b2b');
end


xlabel(['Ratio of betavalues ', new_Strain, '/1a2a']);                              
ylabel('Max peak of infection')
legend('show')

hold off

saveas(gcf, '/Users/u2176312/OneDrive - University of Warwick/Model/New Model/Betatesting/1b2a_infection_200days_ratiobetavalues.pdf')