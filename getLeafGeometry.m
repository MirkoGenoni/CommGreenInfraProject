function leaf = getLeafGeometry(species)
% getLeafGeometry  
%   Return morphological and anatomical parameters
%   for a given plant species (Yamane & Tani, 2024)
%
% Usage:
%   leaf = getLeafGeometry('Quercus acutissima')
%
% Supported species:
%   'Spathiphyllum clevelandii'
%   'Quercus acutissima'
%   'Quercus myrsinaefolia'
%
% Output structure fields:
%   leaf.L_cw        - cell wall thickness (m)
%   leaf.L_ct        - cytosol thickness (m)
%   leaf.L_ias       - intercellular air space path length (m)
%   leaf.Ames_AT     - ratio mesophyll surface area / projected leaf area
%   leaf.f_ias       - fraction of intercellular air space
%   leaf.theta       - fraction of mesophyll adjacent to chloroplasts
%   leaf.Acw_Aias    - ratio of cell wall to intercellular air space area

    switch lower(strtrim(species))
        case {'spathiphyllum clevelandii','spathiphyllum','s. clevelandii'}
            leaf.L_cw      = 3.70e-07;   % m
            leaf.L_ct      = 1.43e-07;   % m
            leaf.L_ias     = 6.38e-05;   % m
            leaf.Ames_AT   = 8.99;       % m2/m2
            leaf.f_ias     = 0.328;      % m3/m3
            leaf.theta     = 0.696;      % m2/m2
            leaf.Acw_Aias  = 0.060;      % m2/m2
        case {'quercus myrsinaefolia','q. myrsinaefolia'}
            leaf.L_cw      = 4.63e-07;   % m
            leaf.L_ct      = 1.44e-07;   % m
            leaf.L_ias     = 8.06e-05;   % m
            leaf.Ames_AT   = 12.22;      % m2/m2
            leaf.f_ias     = 0.319;      % m3/m3      
            leaf.theta     = 0.821;      % m2/m2
            leaf.Acw_Aias  = 0.091;      % m2/m2
        case {'quercus acutissima','q. acutissima'}
            leaf.L_cw      = 3.97e-07;   % m
            leaf.L_ct      = 1.14e-07;   % m
            leaf.L_ias     = 6.36e-05;   % m
            leaf.Ames_AT   = 11.31;      % m2/m2
            leaf.f_ias     = 0.256;      % m3/m3      
            leaf.theta     = 0.820;      % m2/m2
            leaf.Acw_Aias  = 0.117;      % m2/m2
        otherwise
            error('Unknown species: %s\nSupported: S. clevelandii, Q. acutissima, Q. myrsinaefolia', species);
    end
end