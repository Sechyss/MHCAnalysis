proportions = zeros(10, 4);

for i = 1:20
    % Initialize randomNumbers and fourthNumber with zeros to enter the loop
    randomNumbers = zeros(1, 3);
    fourthNumber = 0;
    
    % Generate three non-negative random numbers, ensuring none are zero
    % Generate the fourth number to ensure the sum is 1 and it is non-negative
    while any(randomNumbers == 0) || fourthNumber == 0
        randomNumbers = rand(1, 3);
        fourthNumber = max(0, 1 - sum(randomNumbers));
    end
    
    % Create the vector
    proportions(i, :) = [randomNumbers, fourthNumber];
end

Population = 10000;
betavalues = [0.002 0.002 0.002 0.002];
time_new = 200;
Time_simulation = 400;
new_Strain = '1b2a';
seeds = 123;
parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

peaks = zeros(52,  1);

for i = 1:size(proportions, 1)
    clf;
    ratios = proportions(i, :);
    [out, t] = MultiGenotype_strain_model(Population, ratios, parameters);
    maxpeaks = max(out, [],2);
    peaks = cat(2, peaks, maxpeaks);
end

clf;
peaks(:, 1) = [];

subplot(2,2,1);
hold on;
scatter(proportions(:, 1), peaks(2, :) + peaks(14, :), 'DisplayName','1a2a')
scatter(proportions(:, 1), peaks(3, :) + peaks(15, :), 'DisplayName','1a2b')
scatter(proportions(:, 1), peaks(4, :) + peaks(16, :), 'DisplayName','1b2a')
scatter(proportions(:, 1), peaks(5, :) + peaks(17, :), 'DisplayName','1b2b')
xlabel('Proportions');
ylabel('Max infection');
title('Infection dynamic for Genotype 12');

subplot(2,2,2);
hold on;
scatter(proportions(:, 2), peaks(20, :) + peaks(28, :), 'DisplayName','1a2a')
scatter(proportions(:, 2), peaks(21, :) + peaks(29, :), 'DisplayName','1a2b')
scatter(proportions(:, 2), peaks(22, :) + peaks(30, :), 'DisplayName','1b2a')
scatter(proportions(:, 2), peaks(23, :) + peaks(31, :), 'DisplayName','1b2b')
xlabel('Proportions');
ylabel('Max infection');
title('Infection dynamic for Genotype 11');

subplot(2,2,3);
hold on;
scatter(proportions(:, 3), peaks(34, :) + peaks(42, :), 'DisplayName','1a2a')
scatter(proportions(:, 3), peaks(35, :) + peaks(43, :), 'DisplayName','1a2b')
scatter(proportions(:, 3), peaks(36, :) + peaks(44, :), 'DisplayName','1b2a')
scatter(proportions(:, 3), peaks(37, :) + peaks(45, :), 'DisplayName','1b2b')
xlabel('Proportions');
ylabel('Max infection');
title('Infection dynamic for Genotype 22');

subplot(2,2,4);
hold on;
scatter(proportions(:, 4), peaks(48, :), 'DisplayName','1a2a')
scatter(proportions(:, 4), peaks(49, :), 'DisplayName','1a2b')
scatter(proportions(:, 4), peaks(50, :), 'DisplayName','1b2a')
scatter(proportions(:, 4), peaks(51, :), 'DisplayName','1b2b')
ylabel('Max infection');
title('Infection dynamic for Genotype 33');



