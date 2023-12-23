seeds = 1;

rng(seeds);

proportions = zeros(10, 4);

for i = 1:100
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
Time_simulation = 365*3;
new_Strain = '1a2b';

parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};

peaks = zeros(1, 4);

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

clf;
