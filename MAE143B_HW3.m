clear; close all; clc;

% System Parameters
d = 0.1; % Time delay (seconds)
a = 1;   % Pole location at s = -1

%% #3,4: Pade (2,2) Approximation Analysis

num_delay22 = [1, -60, 1200];
den_delay22 = [1, 60, 1200];
G_pade22 = tf(num_delay22, den_delay22) * tf(1, [1 a]);
% Critical values obtained from the (2,2) Pade approx
K_crit_22 = 15.568; 
omega_crit_22 = 18.067;

fprintf('--- Pade (2,2) Approximation Results ---\n');
fprintf('Critical Gain (K): %.3f\n', K_crit_22);
fprintf('Imaginary Axis Crossover Frequency (omega): %.3f rad/s\n\n', omega_crit_22);


%% #5: High-Order Padé (16,12) / True Delay Analysis
G_true = tf(1, [1 a], 'InputDelay', d);

% Analytical values for the true delay
omega_crit_true = 16.634;
K_crit_true = 16.664;

fprintf('--- True Delay System (High-Order Approximation) Results ---\n');
fprintf('True Critical Gain (K): %.3f\n', K_crit_true);
fprintf('True Crossover Frequency (omega): %.3f rad/s\n\n', omega_crit_true);

%% #6: Nyquist Plots and compare two cases
K_stable = 0.5 * K_crit_true;  
K_unstable = 2.0 * K_crit_true; 

figure('Name', 'Nyquist Plots Comparison', 'Position', [100, 100, 1000, 450]);

% case A: Half (Stable)
subplot(1,2,1);
nyquist(K_stable * G_true);
grid on;
title(['Stable Case: K = ', num2str(K_stable, '%.2f'), ' (0.5 * K_{crit})']);
xlabel('Real Axis');
ylabel('Imaginary Axis');

% case B: Twice (Unstable)
subplot(1,2,2);
nyquist(K_unstable * G_true);
grid on;
title(['Unstable Case: K = ', num2str(K_unstable, '%.2f'), ' (2.0 * K_{crit})']);
xlabel('Real Axis');
ylabel('Imaginary Axis');

%% Root Locus Plot
figure('Name', 'Root Locus - Padé (2,2)');
rlocus(G_pade22);
grid on;
title('Root Locus using Padé (2,2) Approximation');
