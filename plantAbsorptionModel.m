%% VOC Absorption Model based on Yamane & Tani (2024)
clear; clc; close all;

%% Plant and VOC choice
species  = 'spathiphyllum clevelandii';   %'spathiphyllum clevelandii','quercus myrsinaefolia','quercus acutissima'
compound = 'pa';  %'bza','pa','macr','nba','iba','nva','mvk','mek','dek','mnpk','mibk'
%'bza','pa','nba','iba','nva','mek','dek','mnpk','mibk'
%% Environmental conditions (From the experiment)
T  = 298.15;           % temperature (K) [assumed 25 °C]
Pa = 1.013;            % atmospheric pressure (bar) [assumed/approximated]
Vw = 0.890;         % water viscosity @25C [from Hayduk-Laudie]

%% Leaf geometry (from Table 1)
leaf = getLeafGeometry(species);
L_ias = leaf.L_ias;       % Effective diffusion path length in intercellular air spaces (m)
L_cw  = leaf.L_cw;        % Thickness of cell wall (m)
L_ct  = leaf.L_ct;        % Thickness of cytosol (m)
f_ias = leaf.f_ias;       % Fraction of intercellular air space
theta = leaf.theta;       % Fraction of cell wall adjacent to chloroplasts
Ames_AT = leaf.Ames_AT;   % Mesophyll surface area / projected area ratio (m2/m2)
tau = 1.57;               % diffusion path tortuosity [Syvertsen 1995]

%% VOC properties (From supporting materials, Table S1 and S2)
VOC = getVOCproperties(compound);
MB = VOC.MB;             % VOC molecular weight (g/mol) [from standard chemical data]
VprimeB = VOC.VB_prime;  % LeBas volume (cm3/mol)
VB = VOC.VB;        % VOC diffusion volum (cm3/mol)
H  = VOC.H;              % Henry's law constant (mol/kg/bar)
KoW = VOC.Kow;           % octanol-water partition
DA = VOC.DA;             % gas-phase diffusion coefficient (m2/s)
DL = VOC.DL;             % liquid-phase diffusion coefficient (m2/s) 


%% Plant properties (From supporting materials, Table S1 and S2)
plant = getPlantProperties(species, compound);
Ca = plant.Ca;      % Exposure Concentration (ppb)
Ca = Ca*1e-9;    % Convert ppb to mol/m3
A = plant.A;             % Absorption rate (mol/m2/s)
E  = plant.E;            % Transpiration rate (mol/m2/s)
rWb = plant.rWb;         % Leaf boundary layer resistance for water vapor (m2.s/mol)
rWs = plant.rWs;         % Stomatal resistance for water vapor (m2.s/mol)

% Theoritical Values
%Ca = 5e-9;              % Concentration in air 5 ppb (mol/mol) [For experimentation]
%rW_total = 1.2e5;       % Total water vapor resistance [approx from literature‑average]
%rWb = 3e4;              % Leaf boundary layer resistance for water vapor [approx from literature‑average (2-5e4)]
%rWs = rW_total - rWb;   % Stomatal resistance for water vapor [approx from literature‑average (7-8e4)]

%% Diffusion Coefficients
MA = 28.97;    % air molecular weight (g/mol) [Fuller 1966]
VA = 20.1;     % air diffusion volum (cm3/mol) [Fuller 1966]
MW = 18.02;    % water molecular weight (g/mol) [from standard chemical data]
VW = 13.1;     % water diffusion volum (cm3/mol) [from standard chemical data/Reid 1977]

% VOC diffusion in air [Fuller 1966]
% DA = (1e-3 * ((T) .^ 1.75) * sqrt((MA + MB)/(MA*MB))) / ...
%      (Pa * ( (VA^(1/3) + VB^(1/3))^2 ));
% DA = DA * 1e-4;

% Water vapor diffusion [Fuller 1966]
DW = 2.6e-5;  % m2/s for H2O(g) in air @25C
% DW = (1e-3 * ((T) .^ 1.75) * sqrt((MA + MW)/(MA*MW))) / ...
%      (Pa * ( (VA^(1/3) + VW^(1/3))^2 ));
% DW = DW * 1e-4;

% Liquid-phase diffusion (Equation 9) [Hayduk-Laudie]
%DL = 13.26e-9 / (Vw^(1.14) * (VB)^0.589); % (m2/s)
%% Gas-Phase Resistance (Equations S8–S10)

% Equations S8–S9 
rb = rWb * (DW / DA)^(2/3); % Leaf boundary layer resistance [m2.s/mol]
rs = (DW * rWs) / DA; % Stomatal resistance [m2.s/mol]

% ignore ra and rc (text S5-S7)
rg = rb + rs;

rb = rb * 40.9; % Convert from m2.s/mol to s/m2 [multiply by Cair (40 mol/m3)]
rs = rs * 40.9; % Convert from m2.s/mol to s/m2 [multiply by Cair (40 mol/m3)]

% intercellular air space diffusion (Equation 3)
rias = (L_ias * tau) / (DA * f_ias);

R_gas = (rb*40.9) + (rs*40.9) + rias;

%% Concentration of the target VOC partitioned into the liquid phase\

% VOC concentrato=ion at the bottom of the stomata (Equation S5)
CiNum = ((1/rg)-(E/2))*Ca-A;
CiDen = ((1/rg)+(E/2));
Ci = CiNum/CiDen;

% Intercellular concentration (Equation 5)
Cias = Ci - 22.4e-3 * (T) * rias * A / 273.15 ;

% Henry’s law partition (Equation 6 and 7)
PCias = Cias * Pa;           % bar 
Ci_sat = H * PCias * 1e3;    % liquid concentration convertion (mol/m^3)


%% Liquid-Phase Resistance (Equations 8, 10, 11)

% Correction factors (Specified in the paper)
Rf_ct = 0.294;          % Cytosol viscosity/tortuosity correction [Niinemets-Reichstein, 2002]
Rf_cw = 1.0;            % For cell wall (assume no effect)
p_cw  = 0.3;            % Effective porosity of cell wall
p_ct  = 1.0;            % Effective porosity of cytosol (assume no effect)

% Equation 8 Cytosol and Cell Wall
r_cw_prime = L_cw / (Rf_cw * DL * p_cw);
r_ct_prime = L_ct / (Rf_ct * DL * p_ct);

% Equation 10 plasma membrane
r_pl_prime = (VB^0.918) / (6.7e-4 * KoW^0.67);

% Equation 11
AT = 1;  % for normalization
rcw = (AT/(theta*Ames_AT))*r_cw_prime;
rct = (AT/(theta*Ames_AT))*r_ct_prime;
rpl = (AT/(theta*Ames_AT))*r_pl_prime;

% Cytosol concentration (Equation 12)
Cct = Ci_sat - A * (rcw + rpl + rct);
Cmet = 0;  % Concentration in plant [assume fully metabolized]

% Metabolic resistance: dominant term, adjustable (see Fig.6–7)
rmet = (Cct - Cmet) / A;

R_liq = rcw + rpl + rct + rmet;

%% Total Resistance
R_total = R_gas + R_liq;
%R_total_theory = (((Ca-Cias)/(A)) * (273.15/(T*22.4*1e-3)))+((Ci_sat-Cmet)/A);
R_total_theory = ((Ca-Ci)/(A))+(((Ci-Cias)/(A)) * (273.15/(T*22.4*1e-3)))+((Ci_sat-Cmet)/A);


%% Absorption Rate Calculation
%A = (Ca - Cmet) / R_total; % mol/m^2/s^1

%% Concentarion Inside Leaf Calculation

%% ========== OUTPUT SECTION ==========
fprintf('--- Yamane & Tani VOC Absorption Model ---\n');
fprintf('Species: %s | VOC: %s\n', species, compound);
fprintf('Temperature: %.1f K\n', T);
fprintf('Atmospheric concentration (Ca): %.2e mol/mol\n', Ca);
fprintf('\nAbsorption rate, A = %.3e mol m^-2 s^-1\n', A);
fprintf('Cias = %.3e mol/mol\n', Cias);
fprintf('Ci_sat (liquid) = %.3e mol/m^3\n', Ci_sat);
fprintf('Cct   = %.3e mol/m^3\n', Cct);
fprintf('\nResistances (s/m):\n');
fprintf(' rb = %.2e | rs = %.2e | rias = %.2e\n', rb, rs, rias);
fprintf(' rcw = %.2e | rpl = %.2e | rct = %.2e | rmet = %.2e\n', rcw, rpl, rct, rmet);
fprintf(' TOTAL R_total = %.2e\n', R_total);
fprintf(' TOTAL R_total_theory = %.2e\n', R_total_theory);

%% Sensitivity demonstration (Fig.6 behavior)
% rmet_list = logspace(2,6,60);
% A_list = (Ca - Cmet) ./ (R_gas + rcw + rpl + rct + rmet_list);
% figure;
% loglog(rmet_list, A_list, 'b-', 'LineWidth', 1.5);
% xlabel('Metabolic site resistance r_{met} (s/m)');
% ylabel('Absorption rate A (mol m^{-2} s^{-1})');
% title('Effect of r_{met} on VOC absorption (Eq.13)');
% grid on;