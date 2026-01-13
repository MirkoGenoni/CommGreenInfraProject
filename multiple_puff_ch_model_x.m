%clc, clear, close all;
%% Constant system parameters
% emitted VOC mass
       
x_tx = 0;           % tx x-position [m]
y_tx = 0;           % tx y-position [m]
z_tx = 1;           % tx z-position [m]
x_rx = 0.1:0.01:5; % rx x-position [m]
y_rx = [0];    % rx y-position [m] (MISO)
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
    classA_plus = struct('u',0.5,'sigma_y',0.22.*x_rx./(sqrt(1+0.0001.*x_rx)), 'sigma_z',0.2.*x_rx);
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

    M = t/delta_p;

    [C_air_A_plus, C_air_A_far_field_plus] = anisotropic_gaussian_puff(Q, classA_plus, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_A, C_air_A_far_field] = anisotropic_gaussian_puff(Q, classA, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_B, C_air_B_far_field] = anisotropic_gaussian_puff(Q, classB, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_C, C_air_C_far_field] = anisotropic_gaussian_puff(Q, classC, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    [C_air_D, C_air_D_far_field] = anisotropic_gaussian_puff(Q, classD, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_E, C_air_E_far_field] = anisotropic_gaussian_puff(Q, classE, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    %[C_air_F, C_air_F_far_field] = anisotropic_gaussian_puff(Q, classF, t, x_rx, y_rx, z_rx, h,M,delta_p,eps);
    
    % C_air_vector(i) = C_air_D(end);
    % i = i+1;

end

%% Plots 
figure
 plot(x_rx, C_air_A_plus, 'LineWidth', 1.5); hold on
 plot(x_rx, C_air_A_far_field_plus, 'LineWidth', 1.5); 
 plot(x_rx, C_air_D, 'LineWidth', 1.5)
 plot(x_rx, C_air_D_far_field, 'LineWidth', 1.5)
% plot(x_rx, C_air_D, 'LineWidth', 1.5)
% plot(x_rx, C_air_E, 'LineWidth', 1.5)
%plot(x_rx, C_air_F, 'LineWidth', 1.5)
hold off

grid on


% Bigger labels and title (without bold)
xlabel('x [m]', 'FontSize', 15)
ylabel('Concentration [nmol/m^{3}]', 'FontSize', 15)
title('Concentration as function of x', 'FontSize', 16)

% Bigger legend (without bold)
legend({'Class A (0.5 m/s)','Class A (0.5 m/s, far field)','Class D','Class D (far field)'}, 'Location','best', 'FontSize', 14)

% Increase tick labels and axes thickness
set(gca, 'FontSize', 14, 'LineWidth', 1.2)
% figure
%  plot(x_rx, C_air_A_far_field_plus, 'LineWidth', 1.5); hold on
%   plot(x_rx, C_air_A_far_field, 'LineWidth', 1.5); 
% 
%  plot(x_rx, C_air_B_far_field, 'LineWidth', 1.5)
%  plot(x_rx, C_air_C_far_field, 'LineWidth', 1.5)
%  plot(x_rx, C_air_D_far_field, 'LineWidth', 1.5)
% % plot(x_rx, C_air_E, 'LineWidth', 1.5)
% %plot(x_rx, C_air_F, 'LineWidth', 1.5)
% hold off
% 
% grid on
% xlabel('X')
% ylabel('C_{air}')
% legend({'Class A (Far Field)','Class B','Class C','Class D'}, ...
%        'Location','best')
% title('Air Concentration as function of x')

figure
 plot(x_rx, C_air_A_far_field_plus./C_air_A_plus, 'LineWidth', 1.5); hold on
 plot(x_rx, C_air_A_far_field./C_air_A, 'LineWidth', 1.5); 
 plot(x_rx, C_air_B_far_field./C_air_B, 'LineWidth', 1.5)
 plot(x_rx, C_air_C_far_field./C_air_C, 'LineWidth', 1.5)
 plot(x_rx, C_air_D_far_field./C_air_D, 'LineWidth', 1.5)
% plot(x_rx, C_air_E, 'LineWidth', 1.5)
%plot(x_rx, C_air_F, 'LineWidth', 1.5)
hold off
grid on 

% Bigger labels and title (without bold)
xlabel('x [m]', 'FontSize', 15)
ylabel('C_{far field} / C', 'FontSize', 15)
title('C_{far field} / C ratio as function of distance', 'FontSize', 16)

% Bigger legend (without bold)
legend({'Class A (0.5 m/s)','Class A','Class B','Class C','Class D'}, 'Location','best', 'FontSize', 14)

% Increase tick labels and axes thickness
set(gca, 'FontSize', 14, 'LineWidth', 1.2)

%% Peclet
a = 0.06;
b = 0.0015;
c = -1/2;
AA = (a.^2.*x_rx.*((b.*x_rx+1).^(2.*c-1)));
BB = (1+x_rx.*(b.*c+b));
Pe = x_rx./(AA.*BB)
