%clc, clear, close all;
%% Constant system parameters
% emitted VOC mass
       
x_tx = 0;           % tx x-position [m]
y_tx = 0;           % tx y-position [m]
z_tx = 1;           % tx z-position [m]
x_rx = 0.01:0.01:15; % rx x-position [m]
y_rx = [0 0.1 0.2 -0.1 -0.2];    % rx y-position [m] (MISO)
z_rx = 1;           % rx z-position [m]
h = z_tx;          
 
delta_p = 0.0001; 
t_f = 25;
%t = delta_p:delta_p:t_f;
eps = 0;

%Q = 1*ones(1,1000000)*delta_p;  %constant Q
%Q = ((3 - 0.5)*rand(1, 100000))*delta_p;
%Q = linspace(1, 0, 1000000)*delta_p;
Q = integral; % run transmitter_model before

%% BRIGGS stability classes
    classA = struct('u',1,'sigma_y',0.22.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.2.*x_rx);
    classB = struct('u',3,'sigma_y',0.16.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.12.*x_rx);
    classC = struct('u',5,'sigma_y',0.11.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.08.*x_rx./(sqrt(1+0.0002.*x_rx)));
    classD = struct('u',7,'sigma_y',0.08.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.06.*x_rx./(sqrt(1+0.0015.*x_rx)));
    classE = struct('u',10,'sigma_y',0.06.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.03.*x_rx./((1+0.003.*x_rx)));
    classF = struct('u',1,'sigma_y',0.04.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.016.*x_rx./((1+0.003.*x_rx)));


%% Concentration

i = 1;
C_air_vector = 0;
%t_arr = delta_p:delta_p:20;

for t = t_f %do a vector to calculate the concentration of all points at every time

    M = t/delta_p; %number of puffs

    % [C_air_A_SISO, C_air_A_far_field] = anisotropic_gaussian_puff(Q, classA, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % [C_air_B_SISO, C_air_B_far_field] = anisotropic_gaussian_puff(Q, classB, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % [C_air_C_SISO, C_air_C_far_field] = anisotropic_gaussian_puff(Q, classC, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % [C_air_D_SISO, C_air_D_far_field] = anisotropic_gaussian_puff(Q, classD, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);

    [C_air_A_MISO_01m, C_air_A_far_field] = anisotropic_gaussian_puff(Q, classA, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_B_MISO_01m, C_air_B_far_field] = anisotropic_gaussian_puff(Q, classB, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_C_MISO_01m, C_air_C_far_field] = anisotropic_gaussian_puff(Q, classC, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_D_MISO_01m, C_air_D_far_field] = anisotropic_gaussian_puff(Q, classD, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % 
    %[C_air_E, C_air_E_far_field] = anisotropic_gaussian_puff(Q, classE, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_F, C_air_F_far_field] = anisotropic_gaussian_puff(Q, classF, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    
    % C_air_vector(i) = C_air_D(end);
    % i = i+1;

end

%% Plot
figure
h(1) = plot(x_rx, C_air_A_MISO_01m./C_air_A_SISO, 'LineWidth', 1.5); hold on
h(2) = plot(x_rx, C_air_B_MISO_01m./C_air_B_SISO, 'LineWidth', 1.5);
h(3) = plot(x_rx, C_air_C_MISO_01m./C_air_C_SISO, 'LineWidth', 1.5);
h(4) = plot(x_rx, C_air_D_MISO_01m./C_air_D_SISO, 'LineWidth', 1.5);

plot(x_rx, C_air_A_MISO_1m./C_air_A_SISO, '--', 'LineWidth', 1.5, 'Color', get(h(1),'Color'));
plot(x_rx, C_air_B_MISO_1m./C_air_B_SISO, '--', 'LineWidth', 1.5, 'Color', get(h(2),'Color'));
plot(x_rx, C_air_C_MISO_1m./C_air_C_SISO, '--', 'LineWidth', 1.5, 'Color', get(h(3),'Color'));
plot(x_rx, C_air_D_MISO_1m./C_air_D_SISO, '--', 'LineWidth', 1.5, 'Color', get(h(4),'Color'));

hold off
grid on

% Bigger labels and title (without bold)
xlabel('x [m]', 'FontSize', 15)
ylabel('C_{MISO}/C_{SISO}', 'FontSize', 15)
title('Relative MISO BVOCs concentration as function of distance', 'FontSize', 16)

% Bigger legend (without bold)
legend({'Class A','Class B','Class C','Class D'}, 'Location','best', 'FontSize', 14)

% Increase tick labels and axes thickness
set(gca, 'FontSize', 14, 'LineWidth', 1.2)