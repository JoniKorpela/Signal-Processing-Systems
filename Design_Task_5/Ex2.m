clc;
clear;

%  Task 2 

MC = 10^5;
N = 50;
A = 1;
noise_variance = 1;
C_10 = 2;
C_01 = 2;

prior_H0 = 0.5;
prior_H1 = 1 - prior_H0;
S = (ones(1, N) ./ sqrt(N))';

threshold_bayes = A / (2 * sqrt(N));
state_of_nature = NaN(1,MC); 
x_observed_all = zeros(N, MC);
x_observed_means = zeros(1, MC);

PFA_THEORY = qfunc(threshold_bayes / sqrt(noise_variance / N));
PD_THEORY = qfunc((threshold_bayes - A*S(1)) / sqrt(noise_variance / N));
PM_THEORY = 1 - PD_THEORY;
PE_THEORY = PM_THEORY * prior_H1 + PFA_THEORY * prior_H0;
BAYES_RISK_THEORY = C_01 * PM_THEORY * prior_H1 + C_10 * PFA_THEORY * prior_H0;


for Mc = 1:MC

    state_of_nature(Mc) =  binornd(1, prior_H1);

    x =  A * state_of_nature(Mc) .* S + sqrt(noise_variance) * randn(N, 1);

    x_observed_all(:,Mc) = x;

    x_observed_means(Mc) = mean(x);

end




PFA_SIMULATED = mean(x_observed_means(state_of_nature == 0) > threshold_bayes);

PD_SIMULATED = mean(x_observed_means(state_of_nature == 1) > threshold_bayes);

PM_SIMULATED = 1 - PD_SIMULATED;

PE_SIMULATED =  (sum((x_observed_means(state_of_nature == 0) > threshold_bayes) == 1)...
    + sum((x_observed_means(state_of_nature == 1) > threshold_bayes) == 0))...
    / (MC);

BAYES_RISK_SIMULATED = C_01 * PM_SIMULATED * (sum(state_of_nature == 1) / MC)...
    + C_10 * PFA_SIMULATED * (sum(state_of_nature == 0) / MC);

fprintf("PFA_THEORY = %f, PM_THEORY = %f, PE_THEORY = %f\n", PFA_THEORY, PM_THEORY, PE_THEORY);
fprintf("PFA_SIMULATED = %f, PM_SIMULATED = %f, PE_SIMULATED = %f\n", PFA_SIMULATED, PM_SIMULATED, PE_SIMULATED);
fprintf("BAYES_RISK_THEORY = %f\n", BAYES_RISK_THEORY);
fprintf("BAYES_RISK_SIMULATED = %f\n", BAYES_RISK_SIMULATED);


hold on;
h1 = histogram(x_observed_means(state_of_nature == 0), 200, 'Normalization','pdf', 'FaceAlpha', 0.2, 'FaceColor','b', 'EdgeColor','b');
h2 = histogram(x_observed_means(state_of_nature == 1), 200, 'Normalization','pdf', 'FaceAlpha', 0.2, 'FaceColor','r', 'EdgeColor','r');
x =[min(x_observed_means):0.01:max(x_observed_means)];
H0_mean = 0;
H1_mean = A / sqrt(N);
plot(x, normpdf(x, H0_mean, noise_variance / sqrt(N)), 'b', 'LineWidth', 2, 'DisplayName','H_0 theory');
plot(x, normpdf(x, H1_mean, noise_variance / sqrt(N)), 'r', 'LineWidth', 2, 'DisplayName','H_1 theory');
xline(threshold_bayes, 'k', Color='g');
title('PDFs of sample mean x̅ under H_0 and H_1, (MC simulated and theoretical)')
legend('H_0 MC', 'H_1 MC', 'H_0 theory', 'H_1 theory', 'Threshold \gamma');
hold off;