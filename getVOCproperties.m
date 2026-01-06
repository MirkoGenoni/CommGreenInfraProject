function VOC = getVOCproperties(name)
% getVOCproperties  
%   Returns physical constants for VOC from Yamane & Tani (2024)
%   Derived from Supporting Table S2 and literature cited therein.
%
% Usage:
%   VOC = getVOCproperties('Methyl ethyl ketone')
%
% Supported compounds:
%   'Methyl ethyl ketone'
%
% Output structure fields:
%   VOC.MB        - VOC molecular weight (g/mol) [from standard chemical data]
%   VOC.VB_prime  - LeBas volume (cm3/mol)
%   VOC.VB        - VOC diffusion volum (cm3/mol)
%   VOC.DA        - gas-phase diffusion coefficient (m2/s)
%   VOC.H         - Henry's law constant (mol/kg/bar)
%   VOC.DL        - liquid-phase diffusion coefficient (m2/s)
%   VOC.Kow       - octanol-water partition


    switch lower(strtrim(name))
        case {'benzaldehyde','bza'}
            VOC.MB = 106.12;
            VOC.VB_prime = 118; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 8.161e-6;
            VOC.H = 39;
            VOC.DL = 9.109e-10; 
            VOC.Kow = 30.2;
        case {'propionaldehyde','pa'}
            VOC.MB = 58.08;
            VOC.VB_prime = 74; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 1.073e-5;
            VOC.H = 13;
            VOC.DL = 1.2e-9;
            VOC.Kow = 2.02;
        case {'methacrolein','macr'}
            VOC.MB = 70.09;
            VOC.VB_prime = 88.8; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 9.681e-6;
            VOC.H = 6.5;
            VOC.DL = 1.078e-9;
            VOC.Kow = 2.00;
        case {'n-butyraldehyde','nba'}
            VOC.MB = 72.11;
            VOC.VB_prime = 96.2; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 9.331e-6;
            VOC.H = 9.6;
            VOC.DL = 1.028e-9;
            VOC.Kow = 7.59;
        case {'isobutyraldehyde','iba'}
            VOC.MB = 72.11;
            VOC.VB_prime = 96.2; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 9.331e-6;
            VOC.H = 3.3;
            VOC.DL = 1.028e-9;
            VOC.Kow = 6.82;
        case {'n-valeraldehyde','nva'}
            VOC.MB = 86.13;
            VOC.VB_prime = 118; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 8.355e-6;
            VOC.H = 4.4;
            VOC.DL = 9.100e-10;
            VOC.Kow = 23.1;
        case {'methyl vinyl ketone','mvk'}
            VOC.MB = 70.09;
            VOC.VB_prime = 88.8; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 9.681e-6;
            VOC.H = 41;
            VOC.DL = 1.078e-9;
            VOC.Kow = 2.57;
        case {'methyl ethyl ketone','mek'}
            VOC.MB = 72.11;
            VOC.VB_prime = 96.2; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 9.331e-6;
            VOC.H = 20;
            VOC.DL = 1.028e-9;
            VOC.Kow = 1.95;
        case {'diethyl ketone','dek'}
            VOC.MB = 86.13;
            VOC.VB_prime = 118; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 8.355e-6;
            VOC.H = 20;
            VOC.DL = 9.100e-10;
            VOC.Kow = 6.61;
        case {'methyl n-propyl ketone','mnpk'}
            VOC.MB = 86.13;
            VOC.VB_prime = 118; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 8.355e-6;
            VOC.H = 12;
            VOC.DL = 9.100e-10;
            VOC.Kow = 7.08;
        case {'methyl isobutyl ketone','mibk'}
            VOC.MB = 100.16;
            VOC.VB_prime = 141; 
            VOC.VB = 0.875*VOC.VB_prime;
            VOC.DA = 7.626e-6;
            VOC.H = 2.2;
            VOC.DL = 8.224e-10;
            VOC.Kow = 20.4;
        
        otherwise
            error('VOC not listed. Available: Propionaldehyde, n-Butyraldehyde, Methyl ethyl ketone, Benzaldehyde');
    end
end