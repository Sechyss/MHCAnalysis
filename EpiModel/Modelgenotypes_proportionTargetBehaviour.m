strains = {'1a2b'; '1b2a' ;'1b2b'};
for infection= 1:size(strains, 1)
     next_infection = strains{infection};
     disp(['Working on ', next_infection])

    % Preparation of parameters
    
    Combination_genotypes = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];
    collection_of_peaks = cell(1,6);
    collection_of_proportions = cell(1,6);
    
    for i = 1:size(Combination_genotypes, 1)
        combination = Combination_genotypes(i, :);
        seeds = 1; 
        rng(seeds);
        
        % Select the two genotypes to compare
        Genotype1 = combination(1);
        Genotype2 = combination(2);
        
        proportions = zeros (1, 4);
        for y = 1:1000
            % Initialize randomNumbers and create an array with two random numbers
            randomNumbers = zeros(1, 4);
            firstrandomgenotype = rand(1);
            secondrandomgenotype = 1- firstrandomgenotype;
            randomNumbers(Genotype1) = firstrandomgenotype;
            randomNumbers(Genotype2) = secondrandomgenotype;
            
            % Create the vector
            proportions(y, :) = randomNumbers;
        end
        
        Population = 10000;
        betavalues = [0.000015 0.000015 0.000015 0.000015];
        time_new = 200;
        Time_simulation = 365*3;
        new_Strain = next_infection;
        
        parameters = {betavalues, time_new, Time_simulation, new_Strain seeds};
        
        peaks = zeros(1, 4);
        
        % Running of all models with the different proportions
        
        for z = 1:size(proportions, 1)
            clf;
            ratios = proportions(z, :);
            [out, t] = MultiGenotype_strain_model(Population, ratios, {betavalues, time_new, Time_simulation, new_Strain z});
            sum1a2a = max(out(2, :)+ out(14,:)+ out(20, :)+ out(28,:)+ out(34, :)+ out(42,:)+out(48, :));
            sum1a2b = max(out(3, :)+ out(15,:)+ out(21, :)+ out(29,:)+ out(35, :)+ out(43,:)+out(49, :));
            sum1b2a = max(out(4, :)+ out(16,:)+ out(22, :)+ out(30,:)+ out(36, :)+ out(44,:)+out(50, :));
            sum1b2b = max(out(5, :)+ out(17,:)+ out(23, :)+ out(31,:)+ out(37, :)+ out(45,:)+out(51, :));
            maxpeaks = [sum1a2a, sum1a2b, sum1b2a, sum1b2b];
            peaks = cat(1, peaks, maxpeaks);
        end
        
        
        peaks(1, :) = [];
    
        collection_of_proportions {i} = proportions;
        collection_of_peaks {i} = peaks;
    
    end
    
    % Plotting of the increase of the population
    
    figure;
    hold on;
    for i = 1:size(collection_of_peaks,2)
        peaks = collection_of_peaks{i};
        proportions = collection_of_proportions{i};
        combination = Combination_genotypes(i, :);
    
    % Select the column of one of the genotypes to study 1 is Genotype 12, 
    % 2 is 11, 3 is 22, and 4 is 33.
    
        subplot(3,2,i);
        hold on;
        scatter(proportions(:,combination(1)), peaks(:,1), 50, 'red' , 'o', 'DisplayName','1a2a');
        scatter(proportions(:,combination(1)), peaks(:,2), 50, 'green', 'o', 'DisplayName','1a2b');
        scatter(proportions(:,combination(1)), peaks(:,3), 50, 'blue', 'o', 'DisplayName','1b2a');
        scatter(proportions(:,combination(1)), peaks(:,4), 50, 'magenta', 'o', 'DisplayName','1b2b');
    
    
    % Label the xlabel based on the genotypes, 1 is Genotype 12, 2 is 11, 3 is
    % 22, and 4 is 33
        xlabel(['Ratio of Genotype ', num2str(combination(1)), ' & ' , num2str(combination(2))]);
        ylabel('Max peak of infection');
        legend('show');
        hold off
    
    end
    
    
    hold off
    
    saveas(gcf, ['/Users/u2176312/OneDrive - University of Warwick/Model/New Model/Betatesting/',new_Strain, '_infection_200days_Genotype_ratios.pdf'])

end

close('all')