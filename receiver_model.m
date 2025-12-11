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

%% Distance analysis (fig. 5a)
% Tx setup
m_tau_e = 1.1*10^-9;
x_rx = 0:0.01:2;

% Ch and Rx setup
C_L0 = 0;
u_x = 25;             %wind speed x-direction [m/s]
u_y = 0;             %wind speed y-direction [m/s]
u_z = 0;             %wind speed z-direction [m/s]
D = 0.1;             %diffusion coeff. [m^2/s]
x = x_rx;
z = z_rx;
y = y_rx;
h = z_tx;
k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
tau_diffusion = x_rx.^2/D; %time delay in diffusion regime
tau_advection = x_rx./u_x; %time delay in advection regime
tau = (1./tau_diffusion+1./tau_advection).^(-1);  
tau_r = max(tau_advection);

gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
alpha = (1000*P_l*A_l)/(K_lw*M_l);
beta = (alpha*gamma*K_lw)/(1000*K_aw);
delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
C_L = 0.9*exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

mu_noise = -0.1*C_L;          
sigma_noise = abs(mu_noise)/3;
noise = 0;
MontCal = 3000;
for ml = 1:MontCal
    noise = noise + (mu_noise + sigma_noise.*randn(size(C_L))) ;
end
noise = noise/MontCal;
C_L_noise = C_L+noise;
C_L_normalized = C_L_noise./max(C_L_noise);

% Graph
figure; hold on; grid on;
plot(x,C_L_normalized,'LineWidth', 2)
plot(x,0*C_L_normalized+0.55,'LineWidth', 2, 'LineStyle','--')
xlabel('x\_rx');
ylabel('Normalized concentration');
legend('Normalized C_LN', '0.55*MaxC_LN', 'Location', 'best');
title('Distance analysis (fig. 5a)');

%% Distance-delay analysis (fig. 5b)
% Paper setup
x_rx = 0:0.01:2;
u_x = 0.5;
u_x2 = 0.005;
D = 0.1;
Pe = u_x*max(x_rx)/D;
% Graph
tau_diffusion = x_rx.^2/D;
tau_advection = x_rx/u_x;
tau_advection2 = x_rx/u_x2;
tau = (1./tau_diffusion+1./tau_advection2).^(-1);
figure; hold on; grid on;
plot(x_rx, tau_diffusion, 'LineWidth', 2);
plot(x_rx, tau_advection, 'LineWidth', 2);
plot(x_rx, tau, 'LineWidth', 2); %this seems to be wrong, but I think they put the wrong graph in the paper
xlabel('x\_rx [m]');
ylabel('Time delay [s]');
legend('\tau_{diffusion}', '\tau_{advection}', '\tau_{combined}', 'Location', 'best');
title('Distance-delay analysis (fig. 5b)');

%% Distance-Mass analysis (fig. 5c)
% Tx setup
figure; hold on; grid on;
for m_tau_e = [1.1, 3.3, 5.5, 11]*10^-9
    x_rx = 0:0.01:2;
    
    % Ch and Rx setup
    C_L0 = 0;
    u_x = 25;             %wind speed x-direction [m/s]
    u_y = 0;             %wind speed y-direction [m/s]
    u_z = 0;             %wind speed z-direction [m/s]
    D = 0.1;             %diffusion coeff. [m^2/s]
    x = x_rx;
    z = z_rx;
    y = y_rx;
    h = z_tx;
    k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
    tau_diffusion = x_rx.^2./D; %time delay in diffusion regime
    tau_advection = x_rx./u_x; %time delay in advection regime
    tau = (1./tau_diffusion+1./tau_advection).^(-1);  
    tau_r = max(tau_advection);
    gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
    lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
    x = 0:0.01:2;
    alpha = (1000*P_l*A_l)/(K_lw*M_l);
    beta = (alpha*gamma*K_lw)/(1000*K_aw);
    delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
    sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
    C_L = exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

    mu_noise = -0.1*C_L;          
    sigma_noise = abs(mu_noise)/3;
    noise = 0;
    MontCal = 3000;
    for ml = 1:MontCal
        noise = noise + (mu_noise + sigma_noise.*randn(size(C_L))) ;
    end
    noise = noise/MontCal;
    C_L_noise = C_L+noise;
    C_L_normalized = C_L./max(C_L);

    % Graph 
    plot(x,C_L_noise,'LineWidth', 2)
end
xlabel('Distance [m]');
ylabel('Concentration');
legend('m\_tau\_e = 1.1e9', 'm\_tau\_e = 3.3e9', 'm\_tau\_e = 5.5e9','m\_tau\_e = 1.1e8', 'Location', 'best');
title('Distance-Mass analysis (fig. 5c)');

%% Wind speed analysis (fig. 6a)
% Tx setup
figure; hold on; grid on;
for u_x = [1 25 50 100]
    x_rx = 0:0.01:5.5;
    m_tau_e = 1.1*10^-9;
    % Ch and Rx setup
    C_L0 = 0;
    u_y = 0;             %wind speed y-direction [m/s]
    u_z = 0;             %wind speed z-direction [m/s]
    D = 0.1;             %diffusion coeff. [m^2/s]
    x = x_rx;
    z = z_rx;
    y = y_rx;
    h = z_tx;
    k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
    tau_diffusion = x_rx.^2./D; %time delay in diffusion regime
    tau_advection = x_rx./u_x; %time delay in advection regime
    tau = (1./tau_diffusion+1./tau_advection).^(-1);  
    tau_r = max(tau_advection);
    gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
    lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
    alpha = (1000*P_l*A_l)/(K_lw*M_l);
    beta = (alpha*gamma*K_lw)/(1000*K_aw);
    delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
    sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
    C_L = 0.9*exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

    mu_noise = -0.1*C_L;          
    sigma_noise = abs(mu_noise)/3;
    noise = 0;
    MontCal = 3000;
    for ml = 1:MontCal
        noise = noise + (mu_noise + sigma_noise.*randn(size(C_L))) ;
    end
    noise = noise/MontCal;
    C_L_noise = C_L+noise;
    C_L_normalized = C_L_noise./max(C_L_noise);

    % Graph 
    plot(x,C_L_normalized,'LineWidth', 2)
end
ylabel('Normalized concentration');
xlabel('Distance [m]');
legend('u\_x = 1', 'u\_x = 25', 'u\_x = 50','u\_x = 100', 'Location', 'best');
title('Wind speed analysis (fig. 6a)');

%% Wind speed-delay analysis (fig. 6b)
% Paper setup
x_rx = 10;
u_x = logspace(-4,2,20*6); 
D = 0.1;
Pe = u_x*x_rx/D;
% Graph
tau_diffusion = x_rx.^2/D*ones(size(u_x));
tau_advection = x_rx./u_x;
tau = (1./tau_diffusion+1./tau_advection).^(-1);
u_min = min(u_x);
u_max = max(u_x);
idx_diff = u_x >= max(1e-4, u_min) & u_x <= min(1e-3, u_max);   % τ_diffusion
idx_tau  = u_x >  max(1e-3, u_min) & u_x <= min(1e-1, u_max);   % τ
idx_adv  = u_x >= max(1e-1, u_min) & u_x <= min(1e2,  u_max);   % τ_advection
figure; hold on; grid on;
plot(u_x(idx_diff), tau_diffusion(idx_diff), 'LineWidth', 2);
plot(u_x(idx_tau),  tau(idx_tau),          'LineWidth', 2);
plot(u_x(idx_adv),  tau_advection(idx_adv),'LineWidth', 2);
set(gca, 'XScale', 'log');
xlabel('u_x (wind speed) [m/s]');
ylabel('Time delay [s]');
legend('\tau_{diffusion}', '\tau_{combined}', '\tau_{advection}', 'Location', 'best');
title('Wind speed-delay analysis (fig. 6b)');

%% Eddy diffusivity analysis (fig. 6c)
% Tx setup
figure; hold on; grid on;
for D = [0.1 10 35 100]
    u_x = 25;
    x_rx = 0.1:0.01:1.5;
    m_tau_e = 1.1*10^-9;
    % Ch and Rx setup
    C_L0 = 0;
    u_y = 0;             %wind speed y-direction [m/s]
    u_z = 0;             %wind speed z-direction [m/s]
    x = x_rx;
    z = z_rx;
    y = y_rx;
    h = z_tx;
    k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
    tau_diffusion = x_rx.^2./D; %time delay in diffusion regime
    tau_advection = x_rx./u_x; %time delay in advection regime
    tau = (1./tau_diffusion+1./tau_advection).^(-1);  
    tau_r = max(tau_advection);
    gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
    lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
    alpha = (1000*P_l*A_l)/(K_lw*M_l);
    beta = (alpha*gamma*K_lw)/(1000*K_aw);
    delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
    sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
    C_L = 0.9*exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

    mu_noise = -0.1*C_L;          
    sigma_noise = abs(mu_noise)/3;
    noise = 0;
    MontCal = 3000;
    for ml = 1:MontCal
        noise = noise + (mu_noise + sigma_noise.*randn(size(C_L))) ;
    end
    noise = noise/MontCal;
    C_L_noise = C_L+noise;
    C_L_normalized = C_L_noise./max(C_L_noise);

    % Graph 
    plot(x,C_L_normalized,'LineWidth', 2)
    
end
    ylabel('Normalized concentration');
    xlabel('Distance [m]');
    legend('D = 0.1', 'D = 10', 'D = 35','D = 100', 'Location', 'best');
    title('Eddy diffusivity analysis (fig. 6c)');

%% Noise intensity analysis (fig. 6d)
% Tx setup
figure; hold on; grid on;
for nInt = [0 1 3 5 7]
    u_x = 25;
    x_rx = 0.1:0.01:1.5;
    m_tau_e = 1.1*10^-9;
    D = 0.1;
    % Ch and Rx setup
    C_L0 = 0;
    u_y = 0;             %wind speed y-direction [m/s]
    u_z = 0;             %wind speed z-direction [m/s]
    x = x_rx;
    z = z_rx;
    y = y_rx;
    h = z_tx;
    k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
    tau_diffusion = x_rx.^2./D; %time delay in diffusion regime
    tau_advection = x_rx./u_x; %time delay in advection regime
    tau = (1./tau_diffusion+1./tau_advection).^(-1);  
    tau_r = max(tau_advection);
    gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
    lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
    alpha = (1000*P_l*A_l)/(K_lw*M_l);
    beta = (alpha*gamma*K_lw)/(1000*K_aw);
    delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
    sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
    C_L = 0.9*exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

    if nInt==0
        C_L_ref = C_L;
    end

    mu_noise = -0.1*C_L*nInt;          
    sigma_noise = abs(mu_noise)/3;
    noise = 0;
    MontCal = 3000;
    for ml = 1:MontCal
        noise = noise + (mu_noise + sigma_noise.*randn(size(C_L))) ;
    end
    noise = noise/MontCal;
    C_L_noise = C_L+noise;
    C_L_normalized = C_L_noise./max(C_L_ref);

    % Graph 
    plot(x,C_L_normalized,'LineWidth', 2)
    
end
    ylabel('Normalized concentration');
    xlabel('Distance [m]');
    legend('nInt = 0', 'nInt = 0.1', 'nInt = 0.3','nInt = 0.5','nInt = 0.7', 'Location', 'best');
    title('Noise intensity analysis (fig. 6d)');

%% Noise intensity SNR analysis (fig. 6e)
% Tx setup
u_x = 25;
x_rx = 0.75;
m_tau_e = 1.1*10^-9;
D = 0.1;
% Ch and Rx setup
C_L0 = 0;
u_y = 0;             %wind speed y-direction [m/s]
u_z = 0;             %wind speed z-direction [m/s]
x = x_rx;
z = z_rx;
y = y_rx;
h = z_tx;
k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
tau_diffusion = x_rx.^2./D; %time delay in diffusion regime
tau_advection = x_rx./u_x; %time delay in advection regime
tau = (1./tau_diffusion+1./tau_advection).^(-1);  
tau_r = max(tau_advection);
gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
alpha = (1000*P_l*A_l)/(K_lw*M_l);
beta = (alpha*gamma*K_lw)/(1000*K_aw);
delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
C_L = 0.9*exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

noiseIntensity = 1:0.1:10;
SNR = zeros(1,length(noiseIntensity));
nn = 1;

figure; hold on; grid on;
for nInt = 1:0.1:10
    mu_noise = -0.1*C_L*nInt;          
    sigma_noise = abs(mu_noise)/3;

    SNR(1,nn) = 10*log10(((C_L^2)/(sigma_noise^2)));
    nn = nn + 1;
end
    % Graph 
    plot(noiseIntensity,SNR(1,:),'LineWidth', 2)
    ylabel('Normalized concentration');
    xlabel('Distance [m]');
    %legend('nInt = 0', 'nInt = 0.1', 'nInt = 0.3','nInt = 0.5','nInt = 0.7', 'Location', 'best');
    title('Noise intensity SNR analysis (fig. 6e)');

%% Receiver threashold analysis (fig. 6f)
% Tx setup
m_tau_e = 1.1*10^-9;
x_rx = 0:0.01:1;

% Ch and Rx setup
C_L0 = 0;
u_x = 25;             %wind speed x-direction [m/s]
u_y = 0;             %wind speed y-direction [m/s]
u_z = 0;             %wind speed z-direction [m/s]
D = 0.1;             %diffusion coeff. [m^2/s]
x = x_rx;
z = z_rx;
y = y_rx;
h = z_tx;
k = (1/u_x)*D.*x_rx;       %integral of Eddy diffusivities (D = constant)  
tau_diffusion = x_rx.^2/D; %time delay in diffusion regime
tau_advection = x_rx./u_x; %time delay in advection regime
tau = (1./tau_diffusion+1./tau_advection).^(-1);  
tau_r = max(tau_advection);

gamma = (m_tau_e/8).*(1./(pi.*k).^(3/2)); %constant to simplify expressions 
lambda = exp(-(y^2)./(4*k)).*(exp(-((z-h)^2)./(4*k))+exp(-((z+h)^2)./(4*k))); 
alpha = (1000*P_l*A_l)/(K_lw*M_l);
beta = (alpha*gamma*K_lw)/(1000*K_aw);
delta = exp((alpha/(K_lw*M_l*u_x^2)).*(K_lw*M_l*u_x.*x+1000*A_l*P_l*k));
sigma = erf((K_lw*M_l*u_x.*(u_x*tau_r-x)-2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x))+erf((K_lw*M_l*u_x.*x+2000*A_l*P_l*k)./(2*K_lw*M_l*sqrt(k)*u_x));
C_L = 0.9*exp(-alpha*tau_r).*(C_L0+((beta.*lambda.*delta.*sqrt(pi*k))/u_x).*sigma); %0.9

mu_noise = -0.1*C_L;          
sigma_noise = abs(mu_noise)/3;
noise = 0;
MontCal = 3000;
for ml = 1:MontCal
    noise = noise + (mu_noise + sigma_noise.*randn(size(C_L))) ;
end
noise = noise/MontCal;
C_L_noise = C_L+noise;
C_L_normalized = C_L_noise./max(C_L_noise);

% Graph
figure; hold on; grid on;
plot(x,C_L_normalized,'LineWidth', 2)
plot(x,0*C_L_normalized+0.5,'LineWidth', 2, 'LineStyle','--')
plot(x,0*C_L_normalized+0.125,'LineWidth', 2, 'LineStyle','--')
plot(x,0*C_L_normalized+0.05,'LineWidth', 2, 'LineStyle','--')
plot(x,0*C_L_normalized+0.025,'LineWidth', 2, 'LineStyle','--')
xlabel('Distance');
ylabel('Normalized concentration');
legend('Normalized CL', '0.5*maxC', '0.125*maxC','0.05*maxC','0. 025*maxC', 'Location', 'best');
title('Receiver threashold analysis (fig. 6f)');

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