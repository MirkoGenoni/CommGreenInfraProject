%% Plant uptake from air with Gaussian puff channel
%  Wind as carrier + multiple puffs
clear; clc; close all;

%% Time discretization
dt   = 0.1;           % time step [s]
Tend = 600;          % total simulation time [s]
t    = 0:dt:Tend;
Nt   = length(t);

%% Physical plant / leaf parameters (Dummy Values)
A_leaf = 5;        % leaf area [m^2]
V_L    = 0.002;        % effective leaf volume [m^3]
K_L    = 10;           % liquid-side mass transfer coeff [m/s]
a_int  = 1;            % interfacial area density [1/m]
g_s    = 86.4/86400;   % stomatal conductance [m/s]

k_m    = 1/200;       % metabolic rate [1/s]

%% Derived coefficients (FULL EXPRESSIONS)
alpha = (A_leaf * g_s) / (K_L * V_L);   % [1/s]
beta  = (A_leaf * g_s) / V_L;                   % [1/s]
alpha = 1/40; % uptake rate [1/s] 
beta = 0.02; % coupling coefficient 

gamma = alpha + k_m;

%% Receiver state (Using two forms)
C_L_ODE = zeros(1,Nt);    % internal leaf concentration
C_L_SOL = zeros(1,Nt);    % analytical solution

%% Channel / Gaussian puff parameters (Dummy Channel)
U   = 1.0;            % wind speed [m/s]
xL  = 50;             % plant distance from source [m]

M   = 1.0;            % mass per puff
sigma0 = 1.5;         % initial spread [m]
D_turb = 0.05;        % turbulent growth rate

% Puff emission strategy
t0 = 50;              % arrival start time [s]
puff_rate = 0.5;      % puffs per second
Tpuff = 1/puff_rate;

t_emit = t0:Tpuff:(t0+300);
Npuff  = length(t_emit);

% Air concentration at leaf
c_air = zeros(1,Nt);

for n = 1:Nt
    tn = t(n);

    for k = 1:Npuff
        tau = tn - t_emit(k);
        if tau > 0
            sigma = sqrt(sigma0^2 + 2*D_turb*tau);
            x_adv = U*tau;

            c_air(n) = c_air(n) + ...
                M/(sqrt(2*pi)*sigma) * ...
                exp(-(xL - x_adv)^2/(2*sigma^2));
        end
    end
end

% Normalize air concentration
c_air = c_air / max(c_air);

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


