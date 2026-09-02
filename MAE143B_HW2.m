close all; clear; 

% 1. System Settings
d = 12;      % Time delay = 12 seconds (Water takes 12s to travel)
a0 = 0.02;   % System parameter (Speed of water mixing)
g.T = 200;   % Simulation time = 200 seconds

% 2. Define the Plant (The Water Bath)
% We use a math trick called "Pade(2,2)" to model the 12-second time delay.
G = RR_pade(d,2,2) * RR_tf(1, [1/a0 1]); 

% 3. Controller Settings (P Controller)
D = 1;       % Controller gain (K = 1)
P = 1/0.5;   % Prefactor (P = 2) used to fix the final temperature error.

% [Figure 1] Root Locus Plot 
figure(1); 
RR_rlocus(G); 
axis([-.4 .3 -.3 .3]); 
title('Root Locus of G(s)'); 
% Meaning: This shows how the system becomes unstable if the gain K is too high.

% [Figure 2] Bath Temperature y(t)
figure(2); 
RR_step(35 + 10 * P * G * D / (1 + G * D), g); 
axis([0 200 32 55]); 
grid on; 
title('Bath Temperature y(t)'); 
% Meaning: Starts at 35°C and goes up to 45°C. At first look, this graph looks good.

% [Figure 3] Valve Temperature u(t) (Problem!)
figure(3); 
RR_step(35 + 10 * P * D / (1 + G * D), g); 
axis([0 200 40 60]); 
grid on; hold on; 
yline(50, 'r--', 'Max Limit (50°C)'); % Red dashed line at 50°C
title('Valve Temperature u(t)'); 
% Meaning: Look at t = 0. The valve temperature goes up to 55°C! 
% But the maximum possible temperature is 50°C. This is a violation.

%% New design
% Define PI Controller
Kp = 0.6; Ki = 0.012;
D_new = RR_tf([Kp Ki], [1 0]); 
P_new = 1;

%% Test with F_2_2 Pade Approximation
figure(4); RR_step(35 + 10 * P_new * G * D_new / (1 + G * D_new), g); 
axis([0 200 32 55]); grid on; 
title('Figure 4: New Bath Temperature y(t) [Pade 2,2]');

% Figure 5
figure(5); RR_step(35 + 10 * P_new * D_new / (1 + G * D_new), g); 
axis([0 200 32 55]); grid on; hold on; yline(50, 'r--', 'Max Limit (50°C)'); 
title('Figure 5: New Valve Temperature u(t) [Pade 2,2] - SAFE');

%% Test with F_16_13
G16 = RR_pade(d,16,13) * RR_tf(1, [1/a0 1]);

% Figure 6
figure(6); RR_step(35 + 10 * P_new * G16 * D_new / (1 + G16 * D_new), g); 
axis([0 200 32 55]); grid on; 
title('Figure 6: Final Verified Bath Temperature y(t) [Pade 16,13]');

% Figure 7
figure(7); RR_step(35 + 10 * P_new * D_new / (1 + G16 * D_new), g); 
axis([0 200 32 55]); grid on; hold on; yline(50, 'r--', 'Max Limit (50°C)'); 
title('Figure 7: Final Verified Valve Temperature u(t) [Pade 16,13]')