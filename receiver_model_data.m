clc, clear, close all;

% data = load("Rec_05u_5d_5m.mat");
% C_L = data.C_L_ODE;
% C_L_max = max(C_L');
% plot(C_L_max)

dataFiles_1d = ["Rec_05u_1d_5m.mat","Rec_1u_1d_5m.mat","Rec_3u_1d_5m.mat"];
dataFiles_3d = ["Rec_05u_3d_5m.mat","Rec_1u_3d_5m.mat","Rec_3u_3d_5m.mat"];
dataFiles_5d = ["Rec_05u_5d_5m.mat","Rec_1u_5d_5m.mat","Rec_3u_5d_5m.mat"];

td = 0.1:0.1:5;
tdq = linspace(min(td), max(td), 400);


%% Day 1
C_L_1d = zeros(3,length(tdq));
for n = 1:length(dataFiles_1d)
    data = load(dataFiles_1d(n));
    C_L_ODE = data.C_L_ODE;
    C_L_SOL = data.C_L_SOL;
    C_L_SOL_2 = data.C_L_SOL_2;

    C_L_max = max(C_L_SOL');

    C_L_1d(n,:) = interp1(td, C_L_max, tdq, 'pchip');  % shape-preserving
end

%% Day 3
C_L_3d = zeros(3,length(tdq));
for n = 1:length(dataFiles_3d)
    data = load(dataFiles_3d(n));
    C_L_ODE = data.C_L_ODE;
    C_L_SOL = data.C_L_SOL;
    C_L_SOL_2 = data.C_L_SOL_2;

    C_L_max = max(C_L_SOL');

    C_L_3d(n,:) = interp1(td, C_L_max, tdq, 'pchip');  % shape-preserving
end

%% Day 5
C_L_5d = zeros(3,length(tdq));
for n = 1:length(dataFiles_5d)
    data = load(dataFiles_5d(n));
    C_L_ODE = data.C_L_ODE;
    C_L_SOL = data.C_L_SOL;
    C_L_SOL_2 = data.C_L_SOL_2;

    C_L_max = max(C_L_SOL');

    C_L_5d(n,:) = interp1(td, C_L_max, tdq, 'pchip');  % shape-preserving
end

%% Create figure

figure(1); hold on;
% Plot smooth curves
plot(tdq(1:80*2), C_L_1d(1,1:80*2), '-', 'LineWidth', 2);
plot(tdq(1:80*2), C_L_1d(2,1:80*2), '-', 'LineWidth', 2);
plot(tdq(1:80*2), C_L_1d(3,1:80*2), '-', 'LineWidth', 2);
xlim([0 2]);
% Labels and title
xlabel('Distance [m]', 'FontSize', 14)
ylabel('Concentration (C_L) [nmol/m^3]', 'FontSize', 14)
title('Receiver Plant Concentration (Day 1)', 'FontSize', 14)
% Legend
legend({'0.5 m/s','1.0 m/s','3.0 m/s'}, ...
       'Location','northeast','FontSize',14)
% Grid and axes formatting
grid on
set(gca,'FontSize',12)
box on
hold off

figure(2); hold on;
% Plot smooth curves
plot(tdq(1:80*2), C_L_3d(1,1:80*2), '-', 'LineWidth', 2);
plot(tdq(1:80*2), C_L_3d(2,1:80*2), '-', 'LineWidth', 2);
plot(tdq(1:80*2), C_L_3d(3,1:80*2), '-', 'LineWidth', 2);
% % Plot original points with markers
% plot(tdq(1:80*2), C_L_1d(1,1:80*2), '--', 'LineWidth', 1.5)
% plot(tdq(1:80*2), C_L_1d(2,1:80*2), '--', 'LineWidth', 1.5)
% plot(tdq(1:80*2), C_L_1d(3,1:80*2), '--', 'LineWidth', 1.5)
xlim([0 2]);
% Labels and title
xlabel('Distance [m]', 'FontSize', 14)
ylabel('Concentration (C_L) [nmol/m^3]', 'FontSize', 14)
title('Receiver Plant Concentration (Day 3)', 'FontSize', 14)
% Legend
legend({'0.5 m/s','1.0 m/s','3.0 m/s'}, ...
       'Location','northeast','FontSize',14)
% Grid and axes formatting
grid on
set(gca,'FontSize',12)
box on
hold off

figure(3); hold on;
% Plot smooth curves
plot(tdq(1:80*2), C_L_5d(1,1:80*2), '-', 'LineWidth', 2);
plot(tdq(1:80*2), C_L_5d(2,1:80*2), '-', 'LineWidth', 2);
plot(tdq(1:80*2), C_L_5d(3,1:80*2), '-', 'LineWidth', 2);
% % Plot original points with markers
% plot(tdq(1:80*2), C_L_1d(1,1:80*2), '--', 'LineWidth', 1.5)
% plot(tdq(1:80*2), C_L_1d(2,1:80*2), '--', 'LineWidth', 1.5)
% plot(tdq(1:80*2), C_L_1d(3,1:80*2), '--', 'LineWidth', 1.5)
xlim([0 2]);
% Labels and title
xlabel('Distance [m]', 'FontSize', 14)
ylabel('Concentration (C_L) [nmol/m^3]', 'FontSize', 14)
title('Receiver Plant Concentration (Day 5)', 'FontSize', 14)
% Legend
legend({'0.5 m/s','1.0 m/s','3.0 m/s'}, ...
       'Location','northeast','FontSize',14)
% Grid and axes formatting
grid on
set(gca,'FontSize',12)
box on
hold off

%% Data measurments

Diff_at_3d = C_L_3d - C_L_1d;
Diff_at_5d = C_L_5d - C_L_1d;

dd = mean(diff(tdq));
dist = [34 75 115 156];

fprintf("Increase at 0.5 meter for 3d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(1))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_3d(1,round(dist(1))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(1))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_3d(2,round(dist(1))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(1))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_3d(3,round(dist(1))));
fprintf("\n");
fprintf("Increase at 1 meter for 3d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(2))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_3d(1,round(dist(2))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(2))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_3d(2,round(dist(2))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(2))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_3d(3,round(dist(2))));
fprintf("\n");
fprintf("Increase at 1.5 meter for 3d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(3))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_3d(1,round(dist(3))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(3))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_3d(2,round(dist(3))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(3))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_3d(3,round(dist(3))));
fprintf("\n");
fprintf("Increase at 2.0 meter for 3d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(4))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_3d(1,round(dist(4))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(4))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_3d(2,round(dist(4))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(4))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_3d(3,round(dist(4))));

fprintf("Increase at 0.5 meter for 5d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(1))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_5d(1,round(dist(1))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(1))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_5d(2,round(dist(1))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(1))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_5d(3,round(dist(1))));
fprintf("\n");
fprintf("Increase at 1 meter for 5d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(2))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_5d(1,round(dist(2))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(2))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_5d(2,round(dist(2))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(2))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_5d(3,round(dist(2))));
fprintf("\n");
fprintf("Increase at 1.5 meter for 5d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(3))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_5d(1,round(dist(3))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(3))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_5d(2,round(dist(3))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(3))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_5d(3,round(dist(3))));
fprintf("\n");
fprintf("Increase at 2.0 meter for 5d: \n");
fprintf("Wind speed 0.5 m/s, Baseline: %0.3f \n",C_L_1d(1,round(dist(4))));
fprintf("Wind speed 0.5 m/s, diff: %0.3f \n",Diff_at_5d(1,round(dist(4))));
fprintf("Wind speed 1.0 m/s, Baseline: %0.3f \n",C_L_1d(2,round(dist(4))));
fprintf("Wind speed 1.0 m/s, diff: %0.3f \n",Diff_at_5d(2,round(dist(4))));
fprintf("Wind speed 3.0 m/s, Baseline: %0.3f \n",C_L_1d(3,round(dist(4))));
fprintf("Wind speed 3.0 m/s, diff: %0.3f \n",Diff_at_5d(3,round(dist(4))));