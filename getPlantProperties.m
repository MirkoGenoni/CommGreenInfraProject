function plant = getPlantProperties(species, name)
% getVOCproperties  
%   Returns physical constants for plant from Yamane & Tani (2024)
%   Derived from Supporting Table S1 and literature cited therein.
%
% Usage:
%   plant = getPlantProperties(species, VOC)
%
% Supported plants:
%   'Methyl ethyl ketone'
%
% Supported compounds:
%   'Methyl ethyl ketone'
%
% Output structure fields:
%   plant.Ca      - Exposure Concentration (ppb)
%   plant.A       - Absorption rate (mol/m2/s)
%   plant.E       - Transpiration rate (mol/m2/s)
%   plant.rWb     - Leaf boundary layer resistance for water vapor (m2.s/mol)
%   plant.rWs     - Stomatal resistance for water vapor (m2.s/mol)

    switch lower(strtrim(species))
        case {'spathiphyllum clevelandii','spathiphyllum','s. clevelandii'}
             switch lower(strtrim(name))
                 case {'benzaldehyde','bza'}
                    plant.Ca  = 55;
                    plant.A   = 1.9e-10;
                    plant.E   = 2.8e-4;
                    plant.rWb = 2.42;
                    plant.rWs = 44.4;
                case {'propionaldehyde','pa'}
                    plant.Ca  = 25;
                    plant.A   = 3.3e-10;
                    plant.E   = 4.8e-4;
                    plant.rWb = 2.54;
                    plant.rWs = 29.3;
                case {'n-butyraldehyde','nba'}
                    plant.Ca  = 89;
                    plant.A   = 1.3e-9;
                    plant.E   = 6.0e-4;
                    plant.rWb = 2.58;
                    plant.rWs = 21.8;
                case {'isobutyraldehyde','iba'}
                    plant.Ca  = 86;
                    plant.A   = 7.1e-10;
                    plant.E   = 3.6e-4;
                    plant.rWb = 2.53;
                    plant.rWs = 39.4;
                case {'n-valeraldehyde','nva'}
                    plant.Ca  = 97;
                    plant.A   = 8.6e-10;
                    plant.E   = 4.5e-4;
                    plant.rWb = 2.64;
                    plant.rWs = 30.2;
                case {'methyl ethyl ketone','mek'}
                    plant.Ca  = 37;
                    plant.A   = 2.0e-10;
                    plant.E   = 5.4e-4;
                    plant.rWb = 3.23;
                    plant.rWs = 26.5;
                case {'diethyl ketone','dek'}
                    plant.Ca  = 56;
                    plant.A   = 3.3e-10;
                    plant.E   = 4.1e-4;
                    plant.rWb = 2.87;
                    plant.rWs = 25.8;
                case {'methyl n-propyl ketone','mnpk'}
                    plant.Ca  = 60;
                    plant.A   = 3.2e-10;
                    plant.E   = 5.8e-4;
                    plant.rWb = 3.58;
                    plant.rWs = 18.1;
                case {'methyl isobutyl ketone','mibk'}
                    plant.Ca  = 86;
                    plant.A   = 3.6e-10;
                    plant.E   = 4.5e-4;
                    plant.rWb = 2.47;
                    plant.rWs = 16.1;
                otherwise
                    error('VOC not listed. Available: Propionaldehyde, n-Butyraldehyde, Methyl ethyl ketone, Benzaldehyde');
             end
        case {'quercus myrsinaefolia','q. myrsinaefolia'}
            switch lower(strtrim(name))
                case {'methacrolein','macr'}
                    plant.Ca  = 6.4;
                    plant.A   = 1.5e-10;
                    plant.E   = 5.5e-4;
                    plant.rWb = 2.50;
                    plant.rWs = 13.2;
                case {'methyl ethyl ketone','mek'}
                    plant.Ca  = 3.6;
                    plant.A   = 4.2e-11;
                    plant.E   = 5.1e-4;
                    plant.rWb = 2.50;
                    plant.rWs = 21.4;
                otherwise
                    error('VOC not listed. Available: Propionaldehyde, n-Butyraldehyde, Methyl ethyl ketone, Benzaldehyde');
            end
       case {'quercus acutissima','q. acutissima'}
            switch lower(strtrim(name))
                case {'methacrolein','macr'}
                    plant.Ca  = 8.0;
                    plant.A   = 2.0e-10;
                    plant.E   = 1.2e-3;
                    plant.rWb = 2.50;
                    plant.rWs = 10.5;
                case {'methyl ethyl ketone','mek'}
                    plant.Ca  = 3.7;
                    plant.A   = 6.0e-11;
                    plant.E   = 8.4e-4;
                    plant.rWb = 2.50;
                    plant.rWs = 11.8;
                case {'methyl vinyl ketone','mvk'}
                    plant.Ca  = 2.2;
                    plant.A   = 6.8e-11;
                    plant.E   = 6.6e-4;
                    plant.rWb = 0.23;
                    plant.rWs = 91.8;
                case {'mvk_2'}
                    plant.Ca  = 2.8;
                    plant.A   = 4.1e-11;
                    plant.E   = 4.8e-4;
                    plant.rWb = 0.23;
                    plant.rWs = 52.1;
                case {'mvk_3'}
                    plant.Ca  = 3.2;
                    plant.A   = 2.0e-11;
                    plant.E   = 3.1e-4;
                    plant.rWb = 0.23;
                    plant.rWs = 23.1;
                otherwise
                    error('VOC not listed. Available: Propionaldehyde, n-Butyraldehyde, Methyl ethyl ketone, Benzaldehyde');
            end
        otherwise
            error('Unknown species: %s\nSupported: S. clevelandii, Q. acutissima, Q. myrsinaefolia', species);
    end
end