
% Selection of total Population genotype 12

S12_0 = 2600;
I12_1a2a_0 = 10;
I12_1a2b_0 = 1;
I12_1b2a_0 = 1;
I12_1b2b_0 = 2;
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

S11_0 = 2500;
I11_1a2a_0 = 2;
I11_1a2b_0 = 1;
I11_1b2a_0 = 3;
I11_1b2b_0 = 1;
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

S22_0 = 0;
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

S33_0 = 0;
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
    S11_0, I11_1b2a_0, I11_1a2a_0, I11_1a2b_0, I11_1b2b_0, J11_1a2a_0, ...
    J11_1a2b_0, J11_1b2a_0, J11_1b2b_0, R11_1a_0, R11_1b_0, H11_1a_0, ...
    H11_1b_0, M11_0, ...
    ...
    S22_0, I22_1b2a_0, I22_1a2a_0, I22_1a2b_0, I22_1b2b_0,...
    J22_1a2a_0, J22_1a2b_0, J22_1b2a_0, J22_1b2b_0, R22_2a_0, R22_2b_0, ...
    H22_2a_0, H22_2b_0, M22_0,...
    ...
    S33_0, I33_1b2a_0, I33_1a2a_0, I33_1a2b_0, ...
    I33_1b2b_0, M33_0];

% Parameters
gamma_value = 1 / 7;
sigma_value = 1 / 150;
birth_rate = 0.02;
death_rate = 0.00002;
beta_1a2a = 0.0000001;
beta_1a2b = 0.0000003;
beta_1b2a = 0.00000004;
beta_1b2b = 0.00000005;
beta_values = [beta_1a2a, beta_1a2b, beta_1b2a, beta_1b2b];

% Time span for simulation
tspan=[0,800];

% Solve ODE
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t, ret1] = ode45(@(t, y) deriv_equations(t, y, beta_values, gamma_value, sigma_value, death_rate, birth_rate), tspan, y0, options);

function dydt = deriv_equations(t, y, beta_values, gamma, sigma, death, birth)

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
    I11_1b2a = y(20);
    I11_1a2a = y(21);
    I11_1a2b = y(22);
    I11_1b2b = y(23);
    J11_1a2a = y(24);
    J11_1a2b = y(25);
    J11_1b2a = y(26);
    J11_1b2b = y(27);
    R11_1a = y(28);
    R11_1b = y(29);
    H11_1a = y(30);
    H11_1b = y(31);
    M11 = y(32);

    % Genotype 22

    S22 = y(33);
    I22_1b2a = y(34);
    I22_1a2a = y(35);
    I22_1a2b = y(36);
    I22_1b2b = y(37);
    J22_1a2a = y(38);
    J22_1a2b = y(39);
    J22_1b2a = y(40);
    J22_1b2b = y(41);
    R22_2a = y(42);
    R22_2b = y(43);
    H22_2a = y(44);
    H22_2b = y(45);
    M22 = y(46);

    % Genotype 33

    S33 = y(47);
    I33_1b2a = y(48);
    I33_1a2a = y(49);
    I33_1a2b = y(50);
    I33_1b2b = y(51);
    M33 = y(52);

 
  
    % Unpacking of the beta values to estimate lambda force
    beta1a2a = beta_values(1);
    beta1a2b = beta_values(2);
    beta1b2a = beta_values(3);
    beta1b2b = beta_values(4);

    lambda1a2a = beta1a2a * (I12_1a2a + I11_1a2a + I22_1a2a + I33_1a2a + J12_1a2a + J11_1a2a + J22_1a2a);

    lambda1a2b = beta1a2b * (I12_1a2b + I11_1a2b + I22_1a2b + I33_1a2a + J12_1a2b + J11_1a2b + J22_1a2b);

    lambda1b2a = beta1b2a * (I12_1b2a + I11_1b2a + I22_1b2a + I33_1b2a + J12_1b2a + J11_1b2a + J22_1b2a);

    lambda1b2b = beta1b2b * (I12_1b2b + I11_1b2b + I22_1b2b + I33_1b2a + J12_1b2b + J11_1b2b + J22_1b2b);


    tap = sum([S12, I12_1a2a, I12_1b2b, I12_1a2b, I12_1b2a, J12_1a2a, ...
        J12_1a2b, J12_1b2a, J12_1b2b, R12_1a2a, R12_1a2b, R12_1b2a, ...
        R12_1b2b, H12_1a2a, H12_1a2b, H12_1b2a, H12_1b2b, M12, S11, ...
        I11_1a2a, I11_1b2b, I11_1a2b, I11_1b2a, J11_1a2a, J11_1a2b, ...
        J11_1b2a, J11_1b2b, R11_1a, R11_1b, H11_1a, H11_1b, M11, S22, ...
        I22_1a2a, I22_1b2b, I22_1a2b, I22_1b2a, J22_1a2a, J22_1a2b, ...
        J22_1b2a, J22_1b2b, R22_2a, R22_2b, H22_2a, H22_2b, M22, S33, ...
        I33_1a2a, I33_1b2b, I33_1a2b, I33_1b2a, M33]);

    % Equations for genotype 12
    dydt = zeros(52, 1);
    dydt(1) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S12 + birth * tap;
    dydt(2) = -(death + gamma) * I12_1a2a + lambda1a2a * S12;
    dydt(3) = -(death + gamma) * I12_1a2b + lambda1a2b * S12;
    dydt(4) = -(death + gamma) * I12_1b2a + lambda1b2a * S12;
    dydt(5) = -(death + gamma) * I12_1b2b + lambda1b2b * S12;
    dydt(6) = -(death + sigma) * R12_1a2a + gamma * I12_1a2a;
    dydt(7) = -(death + sigma) * R12_1a2b + gamma * I12_1a2b;
    dydt(8) = -(death + sigma) * R12_1b2a + gamma * I12_1b2a;
    dydt(9) = -(death + sigma) * R12_1b2b + gamma * I12_1b2b;
    dydt(10) = -(death + lambda1a2a) * H12_1a2a + sigma * R12_1a2a;
    dydt(11) = -(death + lambda1a2b) * H12_1a2a + sigma * R12_1a2b;
    dydt(12) = -(death + lambda1b2a) * H12_1b2a + sigma * R12_1b2a;
    dydt(13) = -(death + lambda1b2b) * H12_1b2b + sigma * R12_1b2b;
    dydt(14) = -(death + gamma) * J12_1a2a + lambda1a2a * H12_1b2b;
    dydt(15) = -(death + gamma) * J12_1a2b + lambda1a2b * H12_1b2a;
    dydt(16) = -(death + gamma) * J12_1b2a + lambda1b2a * H12_1a2b;
    dydt(17) = -(death + gamma) * J12_1b2b + lambda1b2b * H12_1a2a;
    dydt(18) = -death * M12 + gamma * (J12_1a2a + J12_1b2a + J12_1a2b + J12_1b2b);

    % Equations for genotype 11 (or 13)
    dydt(19) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S11 + birth * tap;
    dydt(20) = -(death + gamma) * I11_1a2a + lambda1a2a * S11;
    dydt(21) = -(death + gamma) * I11_1a2b + lambda1a2b * S11;
    dydt(22) = -(death + gamma) * I11_1b2a + lambda1b2a * S11;
    dydt(23) = -(death + gamma) * I11_1b2b + lambda1b2b * S11;
    dydt(24) = -(death + sigma * R11_1a) + gamma * (I11_1a2a + I11_1a2b);
    dydt(25) = -(death + sigma * R11_1b) + gamma * (I11_1b2a + I11_1b2b);
    dydt(26) = -(death + lambda1a2a + lambda1a2b) * H11_1a + sigma * R11_1a;
    dydt(27) = -(death + lambda1b2a + lambda1b2b) * H11_1a + sigma * R11_1b;
    dydt(28) = -(death + gamma) * J11_1a2a + lambda1a2a * H11_1b;
    dydt(29) = -(death + gamma) * J11_1a2b + lambda1a2b * H11_1b;
    dydt(30) = -(death + gamma) * J11_1b2a + lambda1b2a * H11_1a;
    dydt(31) = -(death + gamma) * J11_1b2b + lambda1b2b * H11_1a;
    dydt(32) = -death * M11 + gamma * (J11_1a2a + J11_1b2a + J11_1a2b + J11_1b2b);

    % Equations for genotype 22 (or 32)
    dydt(33) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S22 + birth * tap;
    dydt(34) = -(death + gamma) * I22_1a2a + lambda1a2a * S22;
    dydt(35) = -(death + gamma) * I22_1a2b + lambda1a2b * S22;
    dydt(36) = -(death + gamma) * I22_1b2a + lambda1b2a * S22;
    dydt(37) = -(death + gamma) * I22_1b2b + lambda1b2b * S22;
    dydt(38) = -(death + sigma * R22_2a) + gamma * (I22_1a2a + I22_1b2a);
    dydt(39) = -(death + sigma * R22_2b) + gamma * (I22_1a2b + I22_1b2b);
    dydt(40) = -(death + lambda1a2a + lambda1b2a) * H22_2a + sigma * R22_2a;
    dydt(41) = -(death + lambda1a2b + lambda1b2b) * H22_2a + sigma * R22_2b;
    dydt(42) = -(death + gamma) * J22_1a2a + lambda1a2a * H22_2b;
    dydt(43) = -(death + gamma) * J22_1a2b + lambda1a2b * H22_2a;
    dydt(44) = -(death + gamma) * J22_1b2a + lambda1b2a * H22_2b;
    dydt(45) = -(death + gamma) * J22_1b2b + lambda1b2b * H22_2a;
    dydt(46) = -death * M22 + gamma * (J22_1a2a + J22_1b2a + J22_1a2b + J22_1b2b);

    % Equations for genotype 33
    dydt(47) = -(death + lambda1a2a + lambda1a2b + lambda1b2a + lambda1b2b) * S33 + birth * 0 + sigma * M33;
    dydt(48) = -(death + gamma) * I33_1a2a + lambda1a2a * S33;
    dydt(49) = -(death + gamma) * I33_1a2b + lambda1a2b * S33;
    dydt(50) = -(death + gamma) * I33_1b2a + lambda1b2a * S33;
    dydt(51) = -(death + gamma) * I33_1b2b + lambda1b2b * S33;
    dydt(52) = -(death + sigma) * M33 + gamma * (I33_1a2a + I33_1b2a + I33_1a2b + I33_1b2b);

end
