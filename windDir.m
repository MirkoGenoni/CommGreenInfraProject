%% Code start
clear; clc; close all;

%% Constant system parameters
M_l = 0.5*10^-3;    %leaf mass [kg] (artemisia tridentata subsp. wyomingensis)
A_l = 25*10^-4;     %leaf area [m^2](artemisia tridentata subsp. wyomingensis)
P_l = 10^-8;        %permeability [m/s]
K_lw = 1;           %partition coeff. between leaf and water
K_aw = 10;          %partition coeff. between air and water
x_tx = 0;           %tx x-position [m]
y_tx = 0;           %tx y-position [m]
z_tx = 1;           %tx z-position [m]
x_rx = 0.001;       %rx x-position [m]
y_rx = 0;           %rx y-position [m]
z_rx = 1;           %rx z-position [m]

%% Model parameters (just a reference to build each section/figure)
% Transmitter (not used)
tau_b = 0;           %time at which constitutive emission rate ends [s]
tau_e = 1;           %time at which constitutive emission rate restarts [s]
t = -2:0.001:2;      %time axis [s]
I = 10*t;            %production/emission rate 
m_tx = (heaviside(t-tau_b)-heaviside(t-tau_e)).*I; %message signal 
m_tau_e = m_tx(t == tau_e);                        %value of the signal at tau_e     
%plot(t,m_tx)
m_tau_e = 1.1*10^-9; %monoterpene emission of Alnus glutinosa under herbivory stress [kg]

% Channel
u_x = 25;            %wind speed x-direction [m/s]
u_y = 0;             %wind speed y-direction [m/s]
u_z = 0;             %wind speed z-direction [m/s]
D = 0.1;             %diffusion coeff. [m^2/s]
x = x_rx;
z = z_rx;
y = y_rx;
Pe = (u_x*x_rx)/D;        %Peclet number, if Pe>>1 => advection regime
k = (1/u_x)*D*x_rx;       %integral of Eddy diffusivities (D = constant)   
h = z_tx;

gamma = (m_tau_e/8)*(1/(pi*k)^(3/2)); %constant to simplify expressions 
lambda = exp(-(y^2)/(4*k))*(exp(-((z-h)^2)/(4*k))+exp(-((z+h)^2)/(4*k))); 
C_air = gamma*exp((-(x-u_x.*t).^2)./(4*k)).*lambda; %concentration C(x,y,z) in air

% Receiver
C_L0 = 0;
tau_diffusion = x_rx.^2/D; %time delay in diffusion regime
tau_advection = x_rx/u_x; %time delay in advection regime
tau = (1./tau_diffusion+1./tau_advection).^(-1); %time delay in mixed regime 
% (however we should use just tau_advection, since the ch model is in the
% advection limite on the x-axis... 
tau_r = max(tau_advection);
alpha = (1000*P_l*A_l)/(K_lw*M_l);
beta = (alpha*gamma*K_lw)/(1000*K_aw);
delta = exp((alpha/(K_lw*M_l*u_x^2))*(K_lw*M_l*u_x*x+1000*A_l*P_l*k));
sigma = erf((K_lw*M_l*u_x*(u_x*tau_r-x)-2000*A_l*P_l*k)/(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x*x+2000*A_l*P_l*k)/(2*K_lw*M_l*sqrt(k)*u_x));
C_L = 0.9*exp(-alpha*tau_r)*(C_L0+beta*(lambda*delta*sqrt(pi*k)/u_x)*sigma);

%% Wind direction analysis (fig. 6d)
x_distances = [0.5 0.8 1.1 1.4 1.7 2.0];   % m
theta_deg = linspace(0,60,100);
theta_rad = deg2rad(theta_deg);

figure(); hold on;

for i = 1:length(x_distances)
    x_r_fixed = x_distances(i);
    tau_r = x_r_fixed / u_x;

    CL = zeros(size(theta_rad));

    for j = 1:length(theta_rad)
        theta = theta_rad(j);

        % 
        x_dash =  x_r_fixed * cos(theta);
        y_dash = -x_r_fixed * sin(theta);

        k = D * x_dash / u_x;

        alpha = (1000*P_l*A_l)/(K_lw*M_l);
        beta = (P_l * A_l * m_tau_e) / (K_aw * M_l * 8 * (pi * k)^(3/2));
        lambda = exp(-((z-h)^2 + y_dash^2) / (4*k)) + ...
                 exp(-((z+h)^2 + y_dash^2) / (4*k));
        Delta = exp((1000 * A_l * P_l * (K_lw * M_l * u_x * x_dash + 1000 * A_l * P_l * k)) / ...
               (K_lw^2 * M_l^2 * u_x^2));

        f1 = ((K_lw * M_l * u_x * (u_x * tau_r - x_dash)) - 2000 * A_l * P_l * k) / ...
             (2 * K_lw * M_l * sqrt(k) * u_x);
        f2 = (K_lw * M_l * u_x * x_dash + 2000 * A_l * P_l * k) / ...
             (2 * K_lw * M_l * sqrt(k) * u_x);

        % Full propagation equation
        CL(j) = exp(-alpha * tau_r) * ...
                ( C_L0 + beta * lambda * Delta * sqrt(pi * k) / u_x * ...
                 (erf(f1) + erf(f2)) );
    end

    % Not considring the noise here
    CLN = CL;

    semilogy(theta_deg, CLN, 'LineWidth',2, 'DisplayName', sprintf('Distance = %.2f m', x_r_fixed));
end

grid on; box on;
xlabel('Degrees of Deviation');
ylabel('C_{L_N} (mol/m^3)');
title('Fig. 7: Wind Direction Analysis');
legend('Location','best');
ylim([1e-100 1]);
set(gca, 'YScale','log', 'LineWidth', 1.2);