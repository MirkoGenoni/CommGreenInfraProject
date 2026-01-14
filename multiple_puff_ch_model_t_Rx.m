clc, clear, close all;
%% Constant system parameters
% emitted VOC mass
       
x_tx = 0;           % tx x-position [m]
y_tx = 0;           % tx y-position [m]
z_tx = 1;           % tx z-position [m]
%x_rx = 0.01:0.01:1; % rx x-position [m]
x_rx = 0.1:0.1:5;%[0.1 1];
y_rx = [0];    % rx y-position [m] (MISO)
z_rx = 1;           % rx z-position [m]
h = z_tx;

dataTx = load("LOX_dr_5d_1e2.mat");
delta_p = dataTx.tau; % 1e-4; 
t_f = dataTx.window_size; % 10;

% delta_p = tau; % 1e-4; 
% t_f = window_size; % 10;

eps = 0; %keep it zero

%Q = [integral, zeros(1,length(integral)/5)]; % run transmitter_model before
Q = dataTx.integral;
%Q = integral;

%% BRIGGS stability classes
    classA = struct('u',1,'sigma_y',0.22.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.2.*x_rx);
    classB = struct('u',3,'sigma_y',0.16.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.12.*x_rx);
    classC = struct('u',5,'sigma_y',0.11.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.08.*x_rx./(sqrt(1+0.0002.*x_rx)));
    classD = struct('u',7,'sigma_y',0.08.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.06.*x_rx./(sqrt(1+0.0015.*x_rx)));
    classE = struct('u',1,'sigma_y',0.06.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.03.*x_rx./((1+0.003.*x_rx)));
    classF = struct('u',1,'sigma_y',0.04.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.016.*x_rx./((1+0.003.*x_rx)));

%% Concentration

i = 1;

t_arr = delta_p:10000000*delta_p:t_f;

C_air_vector_1 = zeros(length(x_rx),length(t_arr));

for t = t_arr %do a vector to calculate the concentration of all points at every time

    M = t/delta_p;

    %[C_air_A, C_air_A_far_field] = anisotropic_gaussian_puff(Q, classA, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_B, C_air_B_far_field] = anisotropic_gaussian_puff(Q, classB, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_C, C_air_C_far_field] = anisotropic_gaussian_puff(Q, classC, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_D, C_air_D_far_field] = anisotropic_gaussian_puff(Q, classD, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_E, C_air_E_far_field] = anisotropic_gaussian_puff(Q, classE, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_F, C_air_F_far_field] = anisotropic_gaussian_puff(Q, classF, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);

    C_air_vector_1(:,i) = C_air_C;

    i = i+1;
    
    t/t_f %completition percentage
end


% Plot for verification
figure
plot(t_arr,C_air_vector_1,'LineWidth', 1.5);
xlabel('t [s]')
ylabel('C_{air} [kg/m^3]')
title('Air Concentration as function of t (x = 1, u = 1)')

figure
plot(Q,'LineWidth', 1.5);
xlabel('t [s]')
ylabel('C_{air} [kg/m^3]')
title('Air Concentration as function of t (x = 1, u = 1)')

nn = length(Q);
nn = nn/3;

avrg = mean(Q(1:nn));
Cplant = avrg*nn;

% save("ClassC_5d_5m_2","C_air_vector_1","t_arr");