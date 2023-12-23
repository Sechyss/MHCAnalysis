function [r_naught] = Contagion_infection (betavalue, population, recovery_rate)
    r_naught = betavalue*population/recovery_rate;
end