%%
close all; clc; clear;

% ---------- LEAF AREA FITTING --------------------------------------------

file = fopen("Data/Transmitter/LOX_to_leaf.csv", "r");
sanizited_data = ...
    fopen("Data/Transmitter/sanitized_LOX_to_leaf_values.csv", "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
end
fclose(file);
fclose(sanizited_data);


data_stressor = ...
    readtable("Data/Transmitter/sanitized_LOX_to_leaf_values.csv", ...
    'PreserveVariableNames',true);

t=data_stressor.x';
func = @(param) param(1)*exp(param(2)*t) - data_stressor.EmissionToLeafArea';

x0=[1e-1,1e-1,1e-1,1e-1,1e-1];
[param_fit,resnorm,residual,eflag,output] = lsqnonlin(func,x0);

original_var = sum((data_stressor.EmissionToLeafArea - ... 
    mean(data_stressor.EmissionToLeafArea)).^2);

r_2_leaf_area = 1-(resnorm/original_var)


t_stressor= linspace(0,30,10000*39);
stressor = param_fit(1)*exp(param_fit(2)*t_stressor);

figure;
scatter(data_stressor.x, data_stressor.EmissionToLeafArea);
hold on; 
plot(t_stressor, stressor);
grid on;
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");
fontsize(16,"points");

%--------------------------------------------------------------------------


% --------- EMISSION FITTING ----------------------------------------------

% DATA SANIFICATION
file = fopen("Data/Transmitter/emission.csv", "r");
sanizited_data = fopen("Data/Transmitter/sanitized_emission_data.csv", ...
    "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
end
fclose(file);
fclose(sanizited_data);

data_emission = readtable("Data/Transmitter/sanitized_emission_data.csv", ...
    'PreserveVariableNames',true);

% EMISSION INTERPOLATION
t=0:5;
emission=data_emission.larvae_2;

dt= 0.001;
t_em_intrp=0:dt:5;

emission_intrp = interp1(t, emission,t_em_intrp,'linear');

% EXTRACT STRESSOR FROM CURRENT EMISSION
leaf_cons = zeros(1,length(t_em_intrp));
for index= 1:length(t_em_intrp)
    temp = find(round(stressor,4)==round(emission_intrp(index),4));
    leaf_cons(1,index) = t_stressor(round((max(temp) + min(temp))/2)-1);
end

% FITTING POLYNOMIAL TO STRESS PROFILE
polynomial_grade = 3:12;
r_2_leaf_area = zeros(1, length(polynomial_grade));
original_var = sum((leaf_cons - ... 
    mean(leaf_cons)).^2);

for k=1:length(polynomial_grade)
    p = polyfit(t_em_intrp, leaf_cons, polynomial_grade(k));
    leaf_fit = polyval(p, t_em_intrp);
    residuals = leaf_cons - leaf_fit;
    resnorm_sq = sum(residuals.^2);
    r_2_leaf_area(k) = 1-(resnorm_sq/original_var);
end

%PLOT STRESSOR PROFILE
figure;
plot(t_em_intrp, leaf_cons);
hold on;
plot(t_em_intrp, leaf_fit);
xlabel("Days");
ylabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
fontsize(16,"points");

%PLOT r^2 VALUE
figure;
plot(polynomial_grade,r_2_leaf_area);
xlabel("Polynomial grade");
ylabel("r^2", 'Interpreter', 'tex');
fontsize(16,"points");

figure;
plot(t_em_intrp,emission_intrp);
hold on;
scatter(t,emission);

% LEAST SQUARE ON DIFFERENTIAL EQUATION
starting_positions= 10.^(-6:-1); %Order of magnitude
%starting_positions=linspace(10^-6,10^-1,1000);
[best_param, error_profile] = ODE_fit(starting_positions,p, ...
    t_em_intrp,emission_intrp);

% SOLVE DIFFERENTIAL EQUATION WITH OPTIMAL PARAMETERS
tsolv=[0 5];
ic = emission_intrp(1);
[t,sol] = ode45(@(t,g) ODE_eq(t,g,best_param(1),best_param(2), ...
    best_param(3), p),tsolv,ic);

% PLOT FITTED SOLUTION TO DIFFERENTIAL EQUATION
plot(t,sol);
xlabel("Days");
ylabel("LOX emission [nmol m^{-2} s^{-1}]", 'Interpreter', 'tex');
fontsize(16,"points");
legend('Emission interpolation', 'Experimental points', ...
    'Emission fitted to differential eq.');

%PLOT ERROR PROFILE OF FITTED DIFF. EQ.
figure;
plot(starting_positions,error_profile);