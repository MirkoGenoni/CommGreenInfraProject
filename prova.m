%% Constant system parameters
clc, clear, close all;

Q = 1;              % emitted VOC mass
x_tx = 0;           % tx x-position [m]
y_tx = 0;           % tx y-position [m]
z_tx = 1;           % tx z-position [m]
x_rx = 0.01:0.01:5; % rx x-position [m]
y_rx = 0;           % rx y-position [m]
z_rx = 1;           % rx z-position [m]
h = z_tx;
t_init = 0.1;

%% Setup Figure and UI
fig = figure('Name', 'Interactive Gaussian Puff', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);
ax = axes('Parent', fig, 'Position', [0.1 0.3 0.8 0.6]);
grid on; hold on;

% Initialize plot lines with dummy data
pA = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class A');
pB = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class B');
pF = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class F');

xlabel('Distance x\_rx [m]');
ylabel('Concentration');
legend('Location', 'best');
title_handle = title(sprintf('Gaussian Puff Distribution at t = %.2f s', t_init));
ylim([0 0.5]); 

%% Add Slider
sld = uicontrol('Parent', fig, 'Style', 'slider', ...
    'Min', 0.1, 'Max', 5, 'Value', t_init, ...
    'Units', 'normalized', 'Position', [0.3 0.1 0.4 0.05], ...
    'Callback', @(src, event) updateSlider(src, pA, pB, pF, title_handle, x_rx, y_rx, z_rx, Q, h));

uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.3 0.16 0.4 0.03], 'String', 'Time (seconds) [t]');

% Initial plot update
updateSlider(sld, pA, pB, pF, title_handle, x_rx, y_rx, z_rx, Q, h);

%% Slider Callback Function
function updateSlider(src, pA, pB, pF, t_title, x_rx, y_rx, z_rx, Q, h)
    t = src.Value; 
    
    % Define stability classes inside callback so they are accessible
    cA.u = 1; 
    cA.sigma_y = 0.22.*x_rx./(sqrt(1+0.0001.*x_rx)); 
    cA.sigma_z = 0.2.*x_rx;
    
    cB.u = 1; 
    cB.sigma_y = 0.16.*x_rx./(sqrt(1+0.0001.*x_rx)); 
    cB.sigma_z = 0.12.*x_rx;
    
    cF.u = 1; 
    cF.sigma_y = 0.04.*x_rx./(sqrt(1+0.0001.*x_rx)); 
    cF.sigma_z = 0.016.*x_rx./((1+0.003.*x_rx));

    % Calculate new concentrations
    C_air_A = anisotropic_gaussian_puff(Q, cA, t, x_rx, y_rx, z_rx, h);
    C_air_B = anisotropic_gaussian_puff(Q, cB, t, x_rx, y_rx, z_rx, h);
    C_air_F = anisotropic_gaussian_puff(Q, cF, t, x_rx, y_rx, z_rx, h);

    % Update plot data
    set(pA, 'YData', C_air_A);
    set(pB, 'YData', C_air_B);
    set(pF, 'YData', C_air_F);
    set(t_title, 'String', sprintf('Gaussian Puff Distribution at t = %.2f s', t));
end
