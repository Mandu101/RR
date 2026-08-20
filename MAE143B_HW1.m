%% MAE143B-HW1
% Yuran_Choi(A18593752)
clear; clc; close all;

%% #1-Lead compensator

D = tf([1 2.68], [1 37.32]); 
figure(1); 
bode(D_lead);
grid on;
title('#1: Bode Plot');

%% #2 Lag and Double Lag

omega_g = 10;          % Wg = 10 rad/s
target_loss_total = -5; % phase: -5 degrees

% 2a: Single Lag
% When z = 100*p, reducing the phase by 5
phase_eq_2a = @(p) atan2d(omega_g, 100*p) - atan2d(omega_g, p) - target_loss_total;

% iteration
p2a = fzero(phase_eq_2a, 0.01); 
z2a = 100 * p2a;

D_lag = tf([1 z2a], [1 p2a]); 

% 2b번 Double Lag
% Each compensator = -2.5 degrees
target_loss_single = target_loss_total / 2; 

% z = 10*p for each
phase_eq_2b = @(p) atan2d(omega_g, 10*p) - atan2d(omega_g, p) - target_loss_single;

% iteration
p2b = fzero(phase_eq_2b, 0.05);
z2b = 10 * p2b;

D_double_lag = tf([1 z2b], [1 p2b])^2; % Square

fprintf('[2a: Single Lag]\n');
fprintf('p = %.4f rad/s\n', p2a);
fprintf('z = %.4f rad/s\n', z2a);
fprintf('--------------------------------------------------\n');
fprintf('[2b: Double Lag]\n');
fprintf('p = %.4f rad/s\n', p2b);
fprintf('z = %.4f rad/s\n', z2b);


% Bode Plot on the same figure
figure(2);

bode(D_lag);          % 2a
hold on;             

bode(D_double_lag);   % 2b on the same plot

grid on;
legend('Single Lag (2a)', 'Double Lag (2b)');
title('#2: Comparison of Lag and Double Lag');

hold off;           

%% #3: Butterworth && Inverse Chebyshev Filter

omega_g = 10;            % Wg = 10 rad/s
target_filter_loss = -5; 

% 3a: Fourth order Butterworth Filter
butter_phase = @(wc) -(atan2d(0.76537*(omega_g/wc), 1 - (omega_g/wc)^2) + ...
    atan2d(1.84776*(omega_g/wc), 1 - (omega_g/wc)^2)) - target_filter_loss;

% Initial Guess = 300
wc_butter = fzero(butter_phase, 300);
D_butter_scaled = tf(wc_butter^4, [1, 2.6131*wc_butter, 3.4142*(wc_butter^2), 2.6131*(wc_butter^3), wc_butter^4]);

% 3b: Inverse Chebyshev filter
cheb2_phase = @(wc) -(atan2d(0.76537*(omega_g/wc), 1.30656 * (1 - (omega_g/wc)^2)) + ...
    atan2d(1.84776*(omega_g/wc), 1.30656 * (1 - (omega_g/wc)^2))) - target_filter_loss;

wc_cheb = fzero(cheb2_phase, 300);

% Calculations
num_cheb = [1, 0, 2.6131*(wc_cheb^2), 0, wc_cheb^4]; 
den_cheb = [1, 2.6131*wc_cheb, 3.4142*(wc_cheb^2), 2.6131*(wc_cheb^3), wc_cheb^4];
D_inv_cheb = tf(num_cheb, den_cheb);

% Results
fprintf('       #3 \n');
fprintf('3a: Butterworth LPF4 = %.4f rad/s\n', wc_butter);
fprintf('3b: Inverse Chebyshev LPF4 = %.4f rad/s\n', wc_cheb);

% Plot
figure(3);
bode(D_butter_scaled); % Draw 3a
hold on;              
bode(D_inv_cheb);      % Draw 3b on the same figure
grid on;
legend('Butterworth (3a)', 'Inverse Chebyshev (3b)');
title('#3: 4th-Order Low-Pass Filters Comparison');
hold off;


%% #4: Loop-Shaping Compensator Assembly and C2D Transformation

% Define the time step and target frequency
H = 0.001;       % Time step (1/1000 seconds)
omega_g = 10;    % Crossover frequency for prewarping (10 rad/s)

% 4a: Assemble the Continuous-Time Compensator
% We combine the Lead, Double Lag, and Butterworth filters together
D_lead = tf([1 2.68], [1 37.32]);             % #1
D_lag_double = tf([1 0.4853], [1 0.04853])^2;   % #2
% D_butter_scaled is already calculated from problem #3

% Total D(s) without gain K
D_s = D_lead * D_lag_double * D_butter_scaled; 

% 4b: Convert D(s) to Digital D(z) using Tustin with Prewarping
opts_c2d = c2dOptions('Method', 'tustin', 'PrewarpFrequency', omega_g);
D_z = c2d(D_s, H, opts_c2d); 

% Extract numerator (b) and denominator (a) coefficients
[num_z, den_z] = tfdata(D_z, 'v');

% Normalize the values so the first term of the denominator is exactly 1
b_raw = num_z / den_z(1);
a_raw = den_z / den_z(1);

% Fix to exactly 4 significant digits
b_coeffs = arrayfun(@(x) str2double(sprintf('%.3e', x)), b_raw);
a_coeffs = arrayfun(@(x) str2double(sprintf('%.3e', x)), a_raw);

% Final results
fprintf('\n#4\n');
fprintf('4a. The total Continuous-Time D(s) has an order of: %d\n\n', order(D_s));
fprintf('4b. Nonzero coefficients rounded to 4 significant digits:\n');

% Print the results
for idx = 2:length(a_coeffs)
    if a_coeffs(idx) ~= 0
        fprintf('  • a_bar_%d (for u[k-%d]) = %0.4g\n', idx-1, idx-1, a_coeffs(idx));
    end
end

% Print numerator coefficients
for idx = 1:length(b_coeffs)
    if b_coeffs(idx) ~= 0
        fprintf('  • b_bar_%d (for e[k-%d]) = %0.4g\n', idx-1, idx-1, b_coeffs(idx));
    end
end

%% #5: Cart-on-a-Hill Plant Simulation

% Define the unstable plant G(s)
numG = 100;
denG = [1 0 -100]; % s^2 - 100
G = tf(numG, denG);

% 5a: Simple Lead Compensator Design
z_simple = 10;
p_simple = 20.4;
D_simple_no_K = tf([1 z_simple], [1 p_simple]);

% Find gain K_simple to satisfy overshoot/rise time at omega = 10
% For simple loop, we evaluate the magnitude at s = j*10
L_simple_no_K = G * D_simple_no_K;
K_simple = 1 / abs(evalfr(L_simple_no_K, 1i*10));
D_simple = K_simple * D_simple_no_K;

% 5b: Loop-Shaping Compensator Gain Tuning
% We use D_s from problem #4 (Lead * Double Lag * Butterworth)
% Tune K_loop to ensure exact crossover at omega_g = 10 rad/s
K_loop = 1 / abs(evalfr(D_s * G, 1i*10));
D_loop_shaping = K_loop * D_s;

% Results
fprintf('#5\n');
fprintf('5a. Simple Lead Design Parameters:\n');
fprintf('  • z = %.2f\n', z_simple);
fprintf('  • p = %.2f\n', p_simple);
fprintf('  • Calculated Gain K_simple = %.4f\n\n', K_simple);
fprintf('5b. Loop-Shaping Tuning Parameters:\n');
fprintf('  • Tuned Overall Gain K_loop = %.4f\n', K_loop);
fprintf('End of Final Report\n\n');

% PLOT 1: ROOT LOCUS
figure(51);
rlocus(G * D_simple);
hold on;
rlocus(G * D_loop_shaping);
grid on;
title('#5: Root Locus Comparison');
legend('Simple Lead System', 'Loop-Shaping System');
hold off;

% PLOT 2: OPEN-LOOP BODE PLOT COMPARISON
figure(52);
bode(G * D_simple);
hold on;
bode(G * D_loop_shaping);
grid on;
title('#5: Open-Loop Bode Plot Comparison');
legend('Simple Lead System', 'Loop-Shaping System');
hold off;

% PLOT 3: CLOSED-LOOP
sys_cl_simple = feedback(G * D_simple, 1);
sys_cl_loop   = feedback(G * D_loop_shaping, 1);

figure(53);
step(sys_cl_simple, 2); % Simple Lead
hold on;
step(sys_cl_loop, 2);   % Loop-Shaping on the same figure
grid on;

ylim([-5, 15]); 

title('#5: Closed-Loop Step Response Comparison');
legend('Simple Lead System (Stable)', 'Loop-Shaping System (Unstable)');
hold off;

