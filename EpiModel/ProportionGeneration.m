function [proportions] = ProportionGeneration(cycles,limit_random, seednumber)

%Creation of proportions for the model and populations, tweak according to
%the number

rng(seednumber);
proportions = zeros(cycles, 4);
for i = 1:cycles
    % Initialize randomNumbers and fourthNumber with zeros to enter the loop
    randomNumbers = zeros(1, 3);
    fourthNumber = 0;
    
    % Generate three non-negative random numbers, ensuring none are zero
    % Generate the fourth number to ensure the sum is 1 and it is non-negative
    while any(randomNumbers == 0) || fourthNumber == 0
        randomNumbers = rand(limit_random, 3);
        fourthNumber = max(0, 1 - sum(randomNumbers));
    end
    
    % Create the vector
    proportions(i, :) = [randomNumbers, fourthNumber];
end
end