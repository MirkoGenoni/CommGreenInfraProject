function new_ch_model()
    %% Constant system parameters
    clc, clear, close all;
    Q = 0.8;              % emitted VOC mass
    x_tx = 0;           % tx x-position [m]
    y_tx = 0;           % tx y-position [m]
    z_tx = 1;           % tx z-position [m]
    x_rx = 0.01:0.01:50; % rx x-position [m]
    y_rx = [0 -0.5 0.5 0.25 -0.25];    % rx y-position [m] (MISO)
    z_rx = 1;           % rx z-position [m]
    h = z_tx;          
    t_init = 1;      

    %% Setup Figure and UI
    fig = figure('Name', 'Interactive Gaussian Puff', ...
                 'Units', 'normalized', ...
                 'Position', [0.2 0.2 0.6 0.6], ...
                 'Color', [1 1 1]); % White background
             
    ax = axes('Parent', fig, 'Position', [0.1 0.3 0.8 0.6]);
    grid on; hold on;

    % Initialize plot lines with zeros
    pA = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class A');
    pB = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class B');
    pC = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class C');
    pD = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class D');
    pE = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class E');
    pF = plot(ax, x_rx, zeros(size(x_rx)), 'LineWidth', 2, 'DisplayName', 'Class F');

    xlabel('Distance x\_rx [m]');
    ylabel('Concentration');
    legend('Location', 'best');
    title_handle = title(sprintf('Gaussian Puff Distribution at t = %.2f s', t_init));
    ylim([0 2]); 

    %% Add Slider 
    sld = uicontrol('Parent', fig, 'Style', 'slider', ...
        'Min', 0.1, 'Max', 60, 'Value', t_init, ...
        'Units', 'normalized', 'Position', [0.3 0.1 0.4 0.05], ...
        'Callback', @(src, event) updateSlider(src, pA, pB, pC, pD, pE, pF, title_handle, x_rx, y_rx, z_rx, Q, h));

    uicontrol('Parent', fig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.3 0.16 0.4 0.03], ...
        'String', 'Move slider to change Time (t)', ...
        'BackgroundColor', [1 1 1]);

    updateSlider(sld, pA, pB, pC, pD, pE, pF, title_handle, x_rx, y_rx, z_rx, Q, h);
end

%% Slider Callback Function
function updateSlider(src, pA, pB, pC, pD, pE, pF, t_title, x_rx, y_rx, z_rx, Q, h)
    t = get(src, 'Value'); % Current slider value
    
    % Recalculate stability classes for the current x_rx
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
    
    % Calculate new concentrations
    C_air_A = anisotropic_gaussian_puff(Q, classA, t, x_rx, y_rx, z_rx, h);
    C_air_B = anisotropic_gaussian_puff(Q, classB, t, x_rx, y_rx, z_rx, h);
    C_air_C = anisotropic_gaussian_puff(Q, classC, t, x_rx, y_rx, z_rx, h);
    C_air_D = anisotropic_gaussian_puff(Q, classD, t, x_rx, y_rx, z_rx, h);
    C_air_E = anisotropic_gaussian_puff(Q, classE, t, x_rx, y_rx, z_rx, h);
    C_air_F = anisotropic_gaussian_puff(Q, classF, t, x_rx, y_rx, z_rx, h);
    
    % Update plot data
    set(pA, 'YData', C_air_A);
    set(pB, 'YData', C_air_B);
    set(pC, 'YData', C_air_C);
    set(pD, 'YData', C_air_D);
    set(pE, 'YData', C_air_E);
    set(pF, 'YData', C_air_F);
    set(t_title, 'String', sprintf('Gaussian Puff Distribution at t = %.2f s', t));
end
