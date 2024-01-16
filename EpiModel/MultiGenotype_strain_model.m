function [outputTable, Timescale] = MultiGenotype_strain_model(population, ratios, params)
    

    % Unpacking of the data 
    
    ratio12 = ratios(1)/sum(ratios);
    ratio11 = ratios(2)/sum(ratios);
    ratio22 = ratios(3)/sum(ratios);
    ratio33 = ratios(4)/sum(ratios);

    beta_values = params{1};
    time_new_strain = params{2};
    time_simulation = params{3};
    new_strain = params{4};
    seed = params{5};

    rng(seed);
    

    % Selection of total Population genotype 12
    
    S12_0 = population * ratio12;
    I12_1a2a_0 = 0;
    I12_1a2b_0 = 0;
    I12_1b2a_0 = 0;
    I12_1b2b_0 = 0;
    J12_1a2a_0 = 0;
    J12_1a2b_0 = 0;
    J12_1b2a_0 = 0;
    J12_1b2b_0 = 0;
    R12_1a2a_0 = 0;
    R12_1a2b_0 = 0;
    R12_1b2a_0 = 0;
    R12_1b2b_0 = 0;
    H12_1a2a_0 = 0;
    H12_1a2b_0 = 0;
    H12_1b2a_0 = 0;
    H12_1b2b_0 = 0;
    M12_0 = 0;
    
    %Selection of total Population genotype 11
    
    S11_0 = population * ratio11;
    I11_1a2a_0 = 0;
    I11_1a2b_0 = 0;
    I11_1b2a_0 = 0;
    I11_1b2b_0 = 0;
    J11_1a2a_0 = 0;
    J11_1a2b_0 = 0;
    J11_1b2a_0 = 0;
    J11_1b2b_0 = 0;
    R11_1a_0 = 0;
    R11_1b_0 = 0;
    H11_1a_0 = 0;
    H11_1b_0 = 0;
    M11_0 = 0;
    
    % Selection of total Population genotype 22
    
    S22_0 = population * ratio22;
    I22_1a2a_0 = 0;
    I22_1a2b_0 = 0;
    I22_1b2a_0 = 0;
    I22_1b2b_0 = 0;
    J22_1a2a_0 = 0;
    J22_1a2b_0 = 0;
    J22_1b2a_0 = 0;
    J22_1b2b_0 = 0;
    R22_2a_0 = 0;
    R22_2b_0 = 0;
    H22_2a_0 = 0;
    H22_2b_0 = 0;
    M22_0 = 0;
    
    % Selection of total Population genotype 33
    
    S33_0 = population * ratio33;
    I33_1a2a_0 = 0;
    I33_1a2b_0 = 0;
    I33_1b2a_0 = 0;
    I33_1b2b_0 = 0;
    M33_0 = 0;
    
    % Initial conditions for genotypes
    y0 = [S12_0, I12_1a2a_0, I12_1a2b_0, I12_1b2a_0, I12_1b2b_0, R12_1a2a_0,... 
        R12_1a2b_0, R12_1b2a_0, R12_1b2b_0, H12_1a2a_0, H12_1a2b_0, H12_1b2a_0,...
        H12_1b2b_0, J12_1a2a_0, J12_1a2b_0, J12_1b2a_0, J12_1b2b_0, M12_0, ...
        ...
        S11_0, I11_1a2a_0, I11_1a2b_0, I11_1b2a_0, I11_1b2b_0, R11_1a_0,...
        R11_1b_0, H11_1a_0, H11_1b_0,...
        J11_1a2a_0, J11_1a2b_0, J11_1b2a_0, J11_1b2b_0, M11_0, ...
        ...
        S22_0, I22_1a2a_0, I22_1a2b_0, I22_1b2a_0, I22_1b2b_0,...
        R22_2a_0, R22_2b_0, H22_2a_0, H22_2b_0, ...
        J22_1a2a_0, J22_1a2b_0, J22_1b2a_0, J22_1b2b_0, M22_0,...
        ...
        S33_0, I33_1a2a_0, I33_1a2b_0, I33_1b2a_0, ...
        I33_1b2b_0, M33_0];
    
    starting_infection = [2 20 34 48];

    randomIndex = randi(length(starting_infection));
    randomNumber = starting_infection(randomIndex);
    y0(randomNumber) = 1; 
    
    % Parameters
    gamma_value = 1 / 14;
    sigma_value = 1 / 180;
    birth_rate = 0.00004;
    death_rate = 0.00004; 
    
    % Time span for simulation
    tspan=0:time_new_strain;
    
    
    % Solve ODE
    options = odeset('AbsTol', 1e-8, 'NonNegative',1:52);
    ode_solution = ode45(@(t, y) deriv_equations(y, beta_values, gamma_value, sigma_value, death_rate, birth_rate), tspan, y0, options);
    
    t = ode_solution.x;
    ret1 = ode_solution.y;

    % Introduction of the new strain

    y0_new = ret1(:, end);

    if strcmpi(new_strain, '1a2b')
        new_infection = [3 21 35 49];

        randomIndex = randi(length(new_infection));
        randomNumber = new_infection(randomIndex);

        y0_new(randomNumber) = 1; 
    
    elseif strcmpi(new_strain, '1b2a')
        new_infection = [4 22 36 50];

        randomIndex = randi(length(new_infection));
        randomNumber = new_infection(randomIndex);

        y0_new(randomNumber) = 1;
    
    elseif strcmpi(new_strain, '1b2b')
        new_infection = [5 23 37 51];

        randomIndex = randi(length(new_infection));
        randomNumber = new_infection(randomIndex);

        y0_new(randomNumber) = 1;
    
    end
    
    % Solve the ODE again using the last entry as the initial condition

    tspan = time_new_strain:time_simulation;
    ode_solution_new = ode45(@(t, y) deriv_equations(y, beta_values, gamma_value, sigma_value, death_rate, birth_rate), tspan, y0_new, options);
    t_new = ode_solution_new.x;
    ret2 = ode_solution_new.y;
    
    % Combine the results from the two simulations
    ret1 = cat(2, ret1, ret2);
    t = cat(2, t, t_new);

    outputTable = ret1;
    Timescale = t;
    
    
    % Plot the results
    legend12 = ["S12", "I12 1a2a ", "I12 1a2b ", "I12 1b2a ", "I12 1b2b ", "R12 1a2a ", ...
        "R12 1a2b ", "R12 1b2a ", "R12 1b2b ", "H12 1a2a ", "H12 1a2b ", "H12 1b2a ",...
        "H12 1b2b ", "J12 1a2a ", "J12 1a2b ", "J12 1b2a ", "J12 1b2b ", "M12 "];
    legend11 = ["S11 ", "I11 1a2a ", "I11 1a2b ", "I11 1b2a ", "I11 1b2b ", "R11 1a ", "R11 1b ",...
        "H11 1a ", "H11 1b ", "J11 1a2a ", "J11 1a2b ", "J11 1b2a ", "J11 1b2b ",...
        "M11 "];
    legend22 = ["S22 ", "I22 1a2a ", "I22 1a2b ", "I22 1b2a ", "I22 1b2b ", "R22 2a ",...
        "R22 2b ", "H22 2a ", "H22 2b ", "J22 1a2a ", "J22 1a2b ", "J22 1b2a ",...
        "J22 1b2b ", "M22 "];
    legend33 = ["S33 ", "I33 1a2a ", "I33 1a2b ", "I33 1b2a ",...
        "I33 1b2b ", "M33 "];
    
    % Create cell array containing all the legend strings
    legends = {legend12, legend11, legend22, legend33};
    
    % Create subplots
    subplot(2, 2, 1);
    hold on;
    plot(t, ret1(1:18, :), 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Population');
    title('Population Dynamics for Genotype 12');
    
    % Add labels to the lines
    legend(legends{1});
    
    subplot(2, 2, 2);
    hold on;
    plot(t, ret1(19:32, :), 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Population');
    title('Population Dynamics for Genotype 11');
    
    % Add labels to the lines
    legend(legends{2});
    
    subplot(2, 2, 3);
    hold on;
    plot(t, ret1(33:46, :), 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Population');
    title('Population Dynamics for Genotype 22');
    
    % Add labels to the lines
    legend(legends{3});
    
    subplot(2, 2, 4);
    hold on;
    plot(t, ret1(47:52, :), 'LineWidth', 2);
    xlabel('Time (days)');
    ylabel('Population');
    title('Population Dynamics for Genotype 33');
    
    % Add labels to the lines
    legend(legends{4});

    hold off;
    
    function dydt = deriv_equations(y, beta_values, gamma, sigma, death, birth)
    
        % Genotype 12
        S12 = y(1);
        I12_1a2a = y(2);
        I12_1a2b = y(3);
        I12_1b2a = y(4);
        I12_1b2b = y(5);
        R12_1a2a = y(6);
        R12_1a2b = y(7);
        R12_1b2a = y(8);
        R12_1b2b = y(9);
        H12_1a2a = y(10);
        H12_1a2b = y(11);
        H12_1b2a = y(12);
        H12_1b2b = y(13);
        J12_1a2a = y(14);
        J12_1a2b = y(15);
        J12_1b2a = y(16);
        J12_1b2b = y(17);
        M12 = y(18);
    
        % Genotype 11
    
        S11 = y(19);
        I11_1a2a = y(20);
        I11_1a2b = y(21);
        I11_1b2a = y(22);
        I11_1b2b = y(23);
        R11_1a = y(24);
        R11_1b = y(25);
        H11_1a = y(26);
        H11_1b = y(27);
        J11_1a2a = y(28);
        J11_1a2b = y(29);
        J11_1b2a = y(30);
        J11_1b2b = y(31);
        M11 = y(32);
    
        % Genotype 22
    
        S22 = y(33);
        I22_1a2a = y(34);
        I22_1a2b = y(35);
        I22_1b2a = y(36);
        I22_1b2b = y(37);
        R22_2a = y(38);
        R22_2b = y(39);
        H22_2a = y(40);
        H22_2b = y(41);
        J22_1a2a = y(42);
        J22_1a2b = y(43);
        J22_1b2a = y(44);
        J22_1b2b = y(45);
        M22 = y(46);
    
        % Genotype 33
    
        S33 = y(47);
        I33_1a2a = y(48);
        I33_1a2b = y(49);
        I33_1b2a = y(50);
        I33_1b2b = y(51);
        M33 = y(52);
    
     
      
        % Unpacking of the beta values to estimate lambda force
        beta1a2a = beta_values(1);
        beta1a2b = beta_values(2);
        beta1b2a = beta_values(3);
        beta1b2b = beta_values(4);
    
        lambda1a2a = beta1a2a * (I12_1a2a + I11_1a2a + I22_1a2a + I33_1a2a + J12_1a2a + J11_1a2a + J22_1a2a);
    
        lambda1a2b = beta1a2b * (I12_1a2b + I11_1a2b + I22_1a2b + I33_1a2b + J12_1a2b + J11_1a2b + J22_1a2b);
    
        lambda1b2a = beta1b2a * (I12_1b2a + I11_1b2a + I22_1b2a + I33_1b2a + J12_1b2a + J11_1b2a + J22_1b2a);
    
        lambda1b2b = beta1b2b * (I12_1b2b + I11_1b2b + I22_1b2b + I33_1b2b + J12_1b2b + J11_1b2b + J22_1b2b);
    
    
        tap = sum(y);
        prob12 = sum(y(1:18))/tap;
        prob11 = sum(y(19:32))/tap;
        prob22 = sum(y(33:46))/tap;
        prob33 = sum(y(37:52))/tap;
    
        % Equations for genotype 12
        dydt = zeros(52, 1);
        dydt(1) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S12 + birth * tap * prob12;
        dydt(2) = -(death + gamma) * I12_1a2a + lambda1a2a * S12;
        dydt(3) = -(death + gamma) * I12_1a2b + lambda1a2b * S12;
        dydt(4) = -(death + gamma) * I12_1b2a + lambda1b2a * S12;
        dydt(5) = -(death + gamma) * I12_1b2b + lambda1b2b * S12;
        dydt(6) = -(death + sigma) * R12_1a2a + gamma * I12_1a2a;
        dydt(7) = -(death + sigma) * R12_1a2b + gamma * I12_1a2b;
        dydt(8) = -(death + sigma) * R12_1b2a + gamma * I12_1b2a;
        dydt(9) = -(death + sigma) * R12_1b2b + gamma * I12_1b2b;
        dydt(10) = -(death + lambda1b2b) * H12_1a2a + sigma * R12_1a2a;
        dydt(11) = -(death + lambda1b2a) * H12_1a2b + sigma * R12_1a2b;
        dydt(12) = -(death + lambda1a2b) * H12_1b2a + sigma * R12_1b2a;
        dydt(13) = -(death + lambda1a2a) * H12_1b2b + sigma * R12_1b2b;
        dydt(14) = -(death + gamma) * J12_1a2a + lambda1a2a * H12_1b2b;
        dydt(15) = -(death + gamma) * J12_1a2b + lambda1a2b * H12_1b2a;
        dydt(16) = -(death + gamma) * J12_1b2a + lambda1b2a * H12_1a2b;
        dydt(17) = -(death + gamma) * J12_1b2b + lambda1b2b * H12_1a2a;
        dydt(18) = -(death * M12) + gamma * (J12_1a2a + J12_1b2a + J12_1a2b + J12_1b2b);
    
        % Equations for genotype 11 (or 13)
        dydt(19) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S11 + birth * tap * prob11;
        dydt(20) = -(death + gamma) * I11_1a2a + lambda1a2a * S11;
        dydt(21) = -(death + gamma) * I11_1a2b + lambda1a2b * S11;
        dydt(22) = -(death + gamma) * I11_1b2a + lambda1b2a * S11;
        dydt(23) = -(death + gamma) * I11_1b2b + lambda1b2b * S11;
        dydt(24) = -(death + sigma * R11_1a) + gamma * (I11_1a2a + I11_1a2b);
        dydt(25) = -(death + sigma * R11_1b) + gamma * (I11_1b2a + I11_1b2b);
        dydt(26) = -(death + lambda1b2a + lambda1b2b) * H11_1a + sigma * R11_1a;
        dydt(27) = -(death + lambda1a2a + lambda1a2b) * H11_1b + sigma * R11_1b;
        dydt(28) = -(death + gamma) * J11_1a2a + lambda1a2a * H11_1b;
        dydt(29) = -(death + gamma) * J11_1a2b + lambda1a2b * H11_1b;
        dydt(30) = -(death + gamma) * J11_1b2a + lambda1b2a * H11_1a;
        dydt(31) = -(death + gamma) * J11_1b2b + lambda1b2b * H11_1a;
        dydt(32) = -(death * M11) + gamma * (J11_1a2a + J11_1b2a + J11_1a2b + J11_1b2b);
    
        % Equations for genotype 22 (or 32)
        dydt(33) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S22 + birth * tap * prob22;
        dydt(34) = -(death + gamma) * I22_1a2a + lambda1a2a * S22;
        dydt(35) = -(death + gamma) * I22_1a2b + lambda1a2b * S22;
        dydt(36) = -(death + gamma) * I22_1b2a + lambda1b2a * S22;
        dydt(37) = -(death + gamma) * I22_1b2b + lambda1b2b * S22;
        dydt(38) = -(death + sigma * R22_2a) + gamma * (I22_1a2a + I22_1b2a);
        dydt(39) = -(death + sigma * R22_2b) + gamma * (I22_1a2b + I22_1b2b);
        dydt(40) = -(death + lambda1a2b + lambda1b2b) * H22_2a + sigma * R22_2a;
        dydt(41) = -(death + lambda1a2a + lambda1b2a) * H22_2a + sigma * R22_2b;
        dydt(42) = -(death + gamma) * J22_1a2a + lambda1a2a * H22_2b;
        dydt(43) = -(death + gamma) * J22_1a2b + lambda1a2b * H22_2a;
        dydt(44) = -(death + gamma) * J22_1b2a + lambda1b2a * H22_2b;
        dydt(45) = -(death + gamma) * J22_1b2b + lambda1b2b * H22_2a;
        dydt(46) = -(death * M22) + gamma * (J22_1a2a + J22_1b2a + J22_1a2b + J22_1b2b);
    
        % Equations for genotype 33
        dydt(47) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S33 + birth * tap * prob33 + sigma * M33;
        dydt(48) = -(death + gamma) * I33_1a2a + lambda1a2a * S33;
        dydt(49) = -(death + gamma) * I33_1a2b + lambda1a2b * S33;
        dydt(50) = -(death + gamma) * I33_1b2a + lambda1b2a * S33;
        dydt(51) = -(death + gamma) * I33_1b2b + lambda1b2b * S33;
        dydt(52) = -(death + sigma) * M33 + gamma * (I33_1a2a + I33_1a2b + I33_1b2a + I33_1b2b);
    
    end
end