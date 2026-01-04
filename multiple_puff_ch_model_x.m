clc, clear, close all;
%% Constant system parameters
% emitted VOC mass
       
x_tx = 0;           % tx x-position [m]
y_tx = 0;           % tx y-position [m]
z_tx = 1;           % tx z-position [m]
x_rx = 0.1:0.01:10; % rx x-position [m]
y_rx = [0];    % rx y-position [m] (MISO)
z_rx = 1;           % rx z-position [m]
h = z_tx;          
 
delta_p = 0.001; 
t_f = 5;
t = delta_p:delta_p:t_f;
eps = 0;

Q = 1*ones(1,1000000)*delta_p;  
%Q = ((3 - 0.5)*rand(1, 100000))*delta_p;
%Q = linspace(1, 0, 1000000)*delta_p;

%% BRIGGS stability classes
    classA = struct('u',1,'sigma_y',0.22.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.2.*x_rx);
    classB = struct('u',1,'sigma_y',0.16.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.12.*x_rx);
    classC = struct('u',1,'sigma_y',0.11.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.08.*x_rx./((1+0.0002.*x_rx)));
    classD = struct('u',1,'sigma_y',0.08.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.06.*x_rx./((1+0.0015.*x_rx)));
    classE = struct('u',1,'sigma_y',0.06.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.03.*x_rx./((1+0.003.*x_rx)));
    classF = struct('u',1,'sigma_y',0.04.*x_rx./(sqrt(1+0.0001.*x_rx)), ...
        'sigma_z',0.016.*x_rx./((1+0.003.*x_rx)));


%% Concentration

i = 1;
C_air_vector = 0;
t_arr = delta_p:10*delta_p:2;

for t = t_f %do a vector to calculate the concentration of all points at every time

    M = t/delta_p;

    %[C_air_A, C_air_A_far_field] = anisotropic_gaussian_puff(Q, classA, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % [C_air_B, C_air_B_far_field] = anisotropic_gaussian_puff(Q, classB, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % [C_air_C, C_air_C_far_field] = anisotropic_gaussian_puff(Q, classC, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_D, C_air_D_far_field] = anisotropic_gaussian_puff(Q, classD, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    % [C_air_E, C_air_E_far_field] = anisotropic_gaussian_puff(Q, classE, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_F, C_air_F_far_field] = anisotropic_gaussian_puff(Q, classF, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    
    % C_air_vector(i) = C_air_D(end);
    % i = i+1;

end


figure
%plot(x_rx, C_air_A, 'LineWidth', 1.5); hold on
% plot(x_rx, C_air_B, 'LineWidth', 1.5)
% plot(x_rx, C_air_C, 'LineWidth', 1.5)
 plot(x_rx, C_air_D, 'LineWidth', 1.5)
% plot(x_rx, C_air_E, 'LineWidth', 1.5)
%plot(x_rx, C_air_F, 'LineWidth', 1.5)
hold off

grid on
xlabel('X')
ylabel('C_{air}')
legend({'Class A','Class B','Class C','Class D','Class E','Class F'}, ...
       'Location','best')
title('Air Concentration as function of x')

% figure
% plot(t_arr,C_air_vector);