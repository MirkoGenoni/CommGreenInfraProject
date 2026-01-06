%% VOC Absorption Model based on Yamane & Tani (2024)
clear; clc; close all;

%% Plant and VOC choice
species  = 'spathiphyllum clevelandii';   %'spathiphyllum clevelandii','quercus myrsinaefolia','quercus acutissima'
compound = ["bza","pa","nba","iba","nva","mek","dek","mnpk","mibk"];  %'bza','pa','macr','nba','iba','nva','mvk','mek','dek','mnpk','mibk'
%'spathiphyllum clevelandii' === ["bza","pa","nba","iba","nva","mek","dek","mnpk","mibk"]
%'quercus myrsinaefolia' === ["macr","mek"]


%% Environmental conditions (From the experiment)
T  = 298.15;           % temperature (K) [assumed 25 °C]
Pa = 1.013;            % atmospheric pressure (bar) [assumed/approximated]
Vw = 0.890;         % water viscosity @25C [from Hayduk-Laudie]
DW = 2.6e-5;  % m2/s for H2O(g) in air @25C
tau = 1.57;               % diffusion path tortuosity [Syvertsen 1995]

% Correction factors (Specified in the paper)
Rf_ct = 0.294;          % Cytosol viscosity/tortuosity correction [Niinemets-Reichstein, 2002]
Rf_cw = 1.0;            % For cell wall (assume no effect)
p_cw  = 0.3;            % Effective porosity of cell wall
p_ct  = 1.0;            % Effective porosity of cytosol (assume no effect)

AT = 1;  % for normalization
Cmet = 0;  % Concentration in plant [assume fully metabolized]

%% Leaf geometry (from Table 1)
leaf = getLeafGeometry(species);
L_ias = leaf.L_ias;       % Effective diffusion path length in intercellular air spaces (m)
L_cw  = leaf.L_cw;        % Thickness of cell wall (m)
L_ct  = leaf.L_ct;        % Thickness of cytosol (m)
f_ias = leaf.f_ias;       % Fraction of intercellular air space
theta = leaf.theta;       % Fraction of cell wall adjacent to chloroplasts
Ames_AT = leaf.Ames_AT;   % Mesophyll surface area / projected area ratio (m2/m2)

rb = zeros(1,length(compound));
rs = zeros(1,length(compound));
rias = zeros(1,length(compound));
rcw = zeros(1,length(compound));
rct = zeros(1,length(compound));
rpl = zeros(1,length(compound));
rmet = zeros(1,length(compound));

for vocNum = 1:length(compound)

    %% VOC properties (From supporting materials, Table S1 and S2)
    VOC = getVOCproperties(compound(vocNum));
    MB = VOC.MB;             % VOC molecular weight (g/mol) [from standard chemical data]
    VprimeB = VOC.VB_prime;  % LeBas volume (cm3/mol)
    VB = VOC.VB;        % VOC diffusion volum (cm3/mol)
    H  = VOC.H;              % Henry's law constant (mol/kg/bar)
    KoW = VOC.Kow;           % octanol-water partition
    DA = VOC.DA;             % gas-phase diffusion coefficient (m2/s)
    DL = VOC.DL;             % liquid-phase diffusion coefficient (m2/s) 
    
    %% Plant properties (From supporting materials, Table S1 and S2)
    plant = getPlantProperties(species, compound(vocNum));
    Ca = plant.Ca;      % Exposure Concentration (ppb)
    Ca = Ca*1e-9;    % Convert ppb to mol/m3
    A = plant.A;             % Absorption rate (mol/m2/s)
    E  = plant.E;            % Transpiration rate (mol/m2/s)
    rWb = plant.rWb;         % Leaf boundary layer resistance for water vapor (m2.s/mol)
    rWs = plant.rWs;         % Stomatal resistance for water vapor (m2.s/mol)
    
    %% Gas-Phase Resistance (Equations S8–S10)
    
    % Equations S8–S9 
    rb(vocNum)= rWb * (DW / DA)^(2/3); % Leaf boundary layer resistance [m2.s/mol]
    rs(vocNum) = (DW * rWs) / DA; % Stomatal resistance [m2.s/mol]
    
    % ignore ra and rc (text S5-S7)
    rg = rb(vocNum) + rs(vocNum);
    rb(vocNum) = rb(vocNum) * 40.9;
    rs(vocNum) = rs(vocNum) * 40.9;
    
    % intercellular air space diffusion (Equation 3)
    rias(vocNum) = (L_ias * tau) / (DA * f_ias);
    
    %% Concentration of the target VOC partitioned into the liquid phase\
    
    % VOC concentrato=ion at the bottom of the stomata (Equation S5)
    CiNum = ((1/rg)-(E/2))*Ca-A;
    CiDen = ((1/rg)+(E/2));
    Ci = CiNum/CiDen;
    
    % Intercellular concentration (Equation 5)
    Cias = Ci - 22.4e-3 * (T) * rias(vocNum) * A / 273.15 ;
    
    % Henry’s law partition (Equation 6 and 7)
    PCias = Cias * Pa;           % bar 
    Ci_sat = H * PCias * 1e3;    % liquid concentration convertion (mol/m^3)
    
    %% Liquid-Phase Resistance (Equations 8, 10, 11)
    
    % Equation 8 Cytosol and Cell Wall
    r_cw_prime = L_cw / (Rf_cw * DL * p_cw);
    r_ct_prime = L_ct / (Rf_ct * DL * p_ct);
    
    % Equation 10 plasma membrane
    r_pl_prime = (VB^0.918) / (6.7e-4 * KoW^0.67);
    
    % Equation 11
    rcw(vocNum) = (AT/(theta*Ames_AT))*r_cw_prime;
    rct(vocNum) = (AT/(theta*Ames_AT))*r_ct_prime;
    rpl(vocNum) = (AT/(theta*Ames_AT))*r_pl_prime;
    
    % Cytosol concentration (Equation 12)
    Cct = Ci_sat - A * (rcw(vocNum) + rpl(vocNum) + rct(vocNum));
    
    % Metabolic resistance: dominant term, adjustable (see Fig.6–7)
    rmet(vocNum) = (Cct - Cmet) / A; 
    rmet(vocNum) = rmet(vocNum);
end

%% Plots settings
figure(1);
hold on;
semilogy(rs,'o-','Color',[0.4 0.8 1],'MarkerFaceColor',[0.4 0.8 1], ...
         'LineWidth',1.2,'MarkerSize',7);       % blue circles (r_s)
semilogy(rb,'s-','Color',[1 0.6 0.6],'MarkerFaceColor',[1 0.6 0.6], ...
         'LineWidth',1.2,'MarkerSize',7);       % pink squares (r_b)
semilogy(rias,'x-','Color',[0.3 1 0.3],'LineWidth',1.2,'MarkerSize',7); % green x (r_ias)

set(gca,'YScale','log','XTick',1:numel(compound),'XTickLabel',compound,...
    'FontSize',12,'LineWidth',1.2);
ylabel('Gas-phase resistance (s m^{-1})');
xlabel('S. clevelandii');
ylim([1e0 1e6]);
legend({'r_s','r_b','r_{ias}'},'Location','northeast');
title('Gas-phase resistances in S. clevelandii');

% section labels (aldehyde/ketone groups)
text(3,2e0,'\leftarrow aldehyde','FontSize',11);
text(7,2e0,'ketone \rightarrow','FontSize',11,'HorizontalAlignment','right');

grid on; 

figure(2);
hold on;
semilogy(rmet,'o-','Color',[1 .5 0],'MarkerFaceColor',[1 .5 0], ...
         'LineWidth',1.2,'MarkerSize',7);                 % orange circles = r_met
semilogy(rpl,'s-','Color',[1 .3 .8],'MarkerFaceColor',[1 .3 .8], ...
         'LineWidth',1.2,'MarkerSize',7);                 % magenta squares = r_pl
semilogy(rct,'^-','Color',[0 0 0],'MarkerFaceColor',[0 0 0], ...
         'LineWidth',1.2,'MarkerSize',7);                 % yellow triangles = r_ct
semilogy(rcw,'^-','Color',[0.8 0.6 1],'MarkerFaceColor',[0.8 0.6 1], ...
         'LineWidth',1.2,'MarkerSize',7);                 % purple triangles = r_cw

set(gca,'YScale','log','XTick',1:numel(compound),'XTickLabel',compound,...
    'FontSize',12,'LineWidth',1.2);
ylabel('Liquid-phase resistance (s m^{-1})');
xlabel('S. clevelandii');
ylim([1e0 1e8]);
legend({'r_{met}','r_{pl}','r_{ct}','r_{cw}'},'Location','northeast');
title('Liquid-phase resistances in S. clevelandii');

% section labels (aldehyde/ketone groups)
text(3,2e0,'\leftarrow aldehyde','FontSize',11);
text(7,2e0,'ketone \rightarrow','FontSize',11,'HorizontalAlignment','right');

grid on; 
