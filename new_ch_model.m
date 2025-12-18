%% Constant system parameters
clc, clear, close all;
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

%% Distance analysis (fig. 5a) !!!NOT READY!!!
% Tx setup
m_tau_e = 1.1*10^-9;
x_rx = 0:0.01:2;

% Ch and Rx setup
C_L0 = 0;

u_x = 25;             % wind speed x-direction [m/s]
u_y = 0;              % wind speed y-direction [m/s]
u_z = 0;              % wind speed z-direction [m/s]
D   = 0.1;            % diffusion coeff. [m^2/s]

x = x_rx;
z = z_rx;
y = y_rx;
h = z_tx;

%Channel anisotropy (y- and z-direction) 
k = (1./u_x).*D.*x_rx;          

%%%
tau_diffusion = x_rx.^2./D;       % time delay in diffusion regime
tau_advection = x_rx./u_x;        % time delay in advection regime
tau = (1./tau_diffusion + 1./tau_advection).^(-1);
tau_r = tau_advection;

gamma = (m_tau_e./8).*(1./(pi.*k).^(3./2));
lambda = exp(-(y.^2)./(4.*k)) .* ...
         (exp(-((z-h).^2)./(4.*k)) + exp(-((z+h).^2)./(4.*k)));
alpha = (1000.*P_l.*A_l)./(K_lw.*M_l);
beta = (alpha.*gamma.*K_lw)./(1000.*K_aw);
delta = exp((alpha./(K_lw.*M_l.*u_x.^2)).* ...
            (K_lw.*M_l.*u_x.*x+1000.*A_l.*P_l.*k));
f1 = ((K_lw.*M_l.*u_x.*(u_x.*tau_r-x))-2000.*A_l.*P_l.*k)./ ...
     (2.*K_lw.*M_l.*sqrt(k).*u_x);
f2 = (K_lw.*M_l.*u_x.*x + 2000.*A_l.*P_l.*k)./ ...
     (2.*K_lw.*M_l.*sqrt(k).*u_x);
sigma = erf(f1) + erf(f2);

C_L = 0.9.*exp(-alpha.*tau_r).* ...
      (C_L0 + ((beta.*lambda.*delta.*sqrt(pi.*k))./u_x).*sigma);

mu_noise = -0.1.*C_L;
sigma_noise = abs(mu_noise)./3;

noise = 0;
MontCal = 3000;

for ml = 1:MontCal
    noise = noise + (mu_noise + sigma_noise.*randn(size(C_L)));
end

noise = noise./MontCal;

C_L_noise = C_L + noise;
C_L_normalized = C_L_noise./max(C_L_noise);

% Graph
figure; hold on; grid on;
plot(x,C_L_normalized,'LineWidth', 2)
plot(x,0*C_L_normalized+0.55,'LineWidth', 2, 'LineStyle','--')
xlabel('x\_rx');
ylabel('Normalized concentration');
legend('Normalized C_LN', '0.55*MaxC_LN', 'Location', 'best');
title('Distance analysis (fig. 5a)');
