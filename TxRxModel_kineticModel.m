%% VOC-based interplant communication: Transmitter + Receiver
%  Implementation from Maitra & Akan (2025)
clear; clc; close all

%% Simulation Parameters
Tend = 200;            % total simulation time [s]
dt = 0.1;              % time step [s]
t = 0:dt:Tend;

%% VOC generation rate due to stress
% transient stress event between 20–60 s
P = zeros(size(t)); 
P(t>=20 & t<=60) = 50;   % [nmol/m2s]

%% Transmitter Parameters (Q. ilex)
ka = 0.002459;          % aqueous pool rate const [1/s]
kl = 3.37791e-5;        % lipid pool rate const [1/s] []
kg = 0.7;               % gas-phase conductance rate [1/s]
n = 0.867;              % partition coeff (water solubility fraction)

% Initial pools
Sa = 0; Sl = 0; Sg = 0;

% results variables
Sa_vec = zeros(size(t));
Sl_vec = zeros(size(t));
Sg_vec = zeros(size(t));
e_tx   = zeros(size(t));   % emission rate [nmol/m2s1]

%% Transmitter Dynamics by integrate ODEs
for i = 2:length(t)
    dSa = n*P(i) - ka*Sa;
    dSl = (1-n)*P(i) - kl*Sl;
    dSg = ka*Sa + kl*Sl - kg*Sg;

    Sa = Sa + dSa*dt;
    Sl = Sl + dSl*dt;
    Sg = Sg + dSg*dt;

    e_tx(i) = kg*Sg;        % emission rate to ambient air
    
    % Pools filling over time
    Sa_vec(i) = Sa;
    Sl_vec(i) = Sl;
    Sg_vec(i) = Sg;
end

%% Link with channel model
% Simulate a simple exponential decay channel as placeholder
channel_gain = exp(-0.01*(1:length(t)));        % dummy attenuation
channel_delay = 20;                             % seconds delay

% Apply delay by shifting emission signal
delay_steps = round(channel_delay/dt);
e_delayed = [zeros(1,delay_steps), e_tx(1:end-delay_steps)];

% Ambient concentration at receiver (example)
c_air = conv(e_delayed, channel_gain, 'same') * dt;
% NOTE: units must be compatible nmol/m^3 etc.

%% Receiver Parameters
Al = 5;                     % leaf area [m2]
Gl = 86.4;                  % leave Conductance []
Gl = Gl / 86400;            % convert from m/day to m/s
Vl = 0.002;                 % leaf volume [m3]
KLA = 10;                   % leaf–air partition coefficient
Pgrowth = 0.035;            % pseudo-first-order rate constant for plant growth
Pgrowth = Pgrowth / 86400;  % convert from 1/day to 1/s

l = (Al*Gl)/(Vl*KLA) + Pgrowth;   % total loss rate [1/s]
gain = (Al*Gl/Vl);              % multiplicative coefficient [1/s]

CR = 0;                          % initial concentration inside leaf [Assume fully metabolized]
CR_vec = zeros(size(t));
b_vec  = zeros(size(t));         % uptake term

%% Receiver Dynamics
for i = 2:length(t)
    b = gain * c_air(i);            % uptake from air
    dCR = -l*CR + b;                % ODE
    CR = CR + dCR*dt;
    CR_vec(i) = CR;
    b_vec(i)  = b;
end

%% OUTPUTS
figure(1)
subplot(3,1,1)
plot(t,P,'k','LineWidth',1.5)
xlabel('Time [s]'); ylabel('P(t) [gen. rate]')
title('Input: VOC generation (stress signal)')

subplot(3,1,2)
plot(t,e_tx,'b','LineWidth',1.5)
xlabel('Time [s]'); ylabel('Emission rate e_{tx}(t)')
title('Transmitter emission to air')

subplot(3,1,3)
plot(t,CR_vec,'r','LineWidth',1.5)
xlabel('Time [s]'); ylabel('C_R(t)')
title('Receiver absorbed concentration in leaf')

sgtitle('VOC-based Interplant MC: Transmitter–Channel–Receiver')

%% Decision Logic (Detecting stress)
% Threshold-based detection
threshold = 0.5 * max(CR_vec);
is_stress_detected = CR_vec > threshold;