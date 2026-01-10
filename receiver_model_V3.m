%% VOC Absorption Model based on Yamane & Tani (2024)
clear; clc; close all;

%% Channel parameters (from 'multiple_puff_ch_model_t.m')
c_air_2sec = load('Data/Channel/concentration_delta_t_0.001s_2_seconds_x_1_m.mat');
c_air_5days = load('Data/Channel/concentration_delta_t_1000s_5_days_x_01_m.mat');
c_air_5days_1m = load('Data/Channel/concentration_delta_t_1000s_5_days_x_1_m.mat');
c_air = c_air_5days_1m.C_air_vector_1;
t_2Sec = load('Data/Channel/time_delta_t_0.001s_2_seconds_x_1_m.mat');
t_5Day = load('Data/Channel/time_delta_t_1000s_5_days.mat');
%t = t_2Sec.t_arr;
t = t_5Day.t_arr;
dt = mean(diff(t));
Nt   = length(t);

% Plot for verification
figure
plot(t,c_air,'LineWidth', 1.5);
xlabel('t [s]')
ylabel('C_{air} [kg/m^3]')
title('Air Concentration as function of t (x = 1, u = 1)')

%% Environmental conditions
T = 301.15;      % K (28C)
P = 101325;      % Pa
R = 8.314;       % J/mol/K

%% Physical plant / leaf parameters

A_leaf = 5e-4;               % leaf area [m^2] (3 – 7 × 10⁻⁴ m² range) (https://gobotany.nativeplanttrust.org/species/alnus/glutinosa/)
V_L    = 2e-3;               % effective leaf volume [m^3]
K_LA    = 500;               % liquid-side mass transfer coeff [m/s]
%K_L    = 1;                 % For experimentation;
g_s_mol    = 0.85;           % stomatal conductance [m/s]
g_s = g_s_mol * R * T / P;   % convertion to m/s
k_m    = 0;              % metabolic rate [1/s]

%% Derived coefficients (FULL EXPRESSIONS)
alpha = (A_leaf * g_s) / (K_LA * V_L);   % [1/s]
beta  = (A_leaf * g_s) / V_L;            % [1/s]
gamma = alpha + k_m;

%% Receiver state
C_L_ODE = zeros(1,Nt);
C_L_SOL = zeros(1,Nt);

%% Plant uptake ODE (implicit Euler)
for n = 1:Nt-1
    C_L_ODE(n+1) = ...
        ( C_L_ODE(n) + dt*beta*c_air(n+1) ) / ...
        ( 1 + dt*gamma );
end

%% Plant uptake analytical solution (convolution)
for n = 2:Nt
    kernel = exp(-gamma*(t(n) - t(1:n)));
    C_L_SOL(n) = beta * trapz(t(1:n), kernel .* c_air(1:n));
end

%% Saturation Condition and Time to saturation
% Saturation
C_L_SS = (beta/gamma) * max(c_air);
fprintf("For constant air concentration A0 = %f \n", max(c_air));
fprintf("The saturation concentration in the leaf CL_SS = %f \n", C_L_SS);
% Time to certain leaf concentration
C_L_targ =  C_L_SS*0.5;
t_sat = -(1/gamma)*log(1-(C_L_targ/C_L_SS));
fprintf("Time to reach 50%% of internal concentration: \n");
fprintf("%d hours \n", t_sat/(60*60));

%% Ploting the results
figure; hold on; grid on;

subplot(1,2,1);

yyaxis left
plot(t, c_air, 'LineWidth', 2)
ylabel('Air concentration at leaf')

yyaxis right
plot(t, C_L_ODE, 'LineWidth', 2)
ylabel('Leaf internal concentration')

xlabel('Time [s]')
title('Gaussian Puff Channel + Plant Uptake (ODE)')

legend('Air concentration', 'Leaf concentration', 'Location', 'southeast')

subplot(1,2,2);

yyaxis left
plot(t, c_air, 'LineWidth', 2)
ylabel('Air concentration at leaf')

yyaxis right
plot(t, C_L_SOL, 'LineWidth', 2)
ylabel('Leaf internal concentration')

xlabel('Time [s]')
title('Gaussian Puff Channel + Plant Uptake (Analytical)')

legend('Air concentration', 'Leaf concentration', 'Location', 'southeast')
