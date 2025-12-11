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
func = @(param) (param(1)+param(2).*t+param(3).*t.^2+param(4).*t.^3 + ... 
    param(5).*t.^4) - data_stressor.EmissionToLeafArea';

x0=[1e-1,1e-1,1e-1,1e-1,1e-1];
[param_fit,resnorm,residual,eflag,output] = lsqnonlin(func,x0);

original_var = sum((data_stressor.EmissionToLeafArea - ... 
    mean(data_stressor.EmissionToLeafArea)).^2);

r_2_leaf_area = 1-(resnorm/original_var)


t_stressor= linspace(0,30,10000*39);
stressor = param_fit(1)+param_fit(2).*t_stressor+param_fit(3).*t_stressor.^2+ ...
    param_fit(4).*t_stressor.^3 + param_fit(5).*t_stressor.^4;

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

% INTERPOLATION
t=0:5;
emission=data_emission.larvae_2;

dt= 1;
t_em_intrp=0:dt:5;

emission_intrp = interp1(t, emission,t_em_intrp,'linear');

figure;
plot(t_em_intrp,emission_intrp);
hold on;
scatter(t,emission);

x0=[1,1,1];
[param_fit,resnorm,residual,eflag,output] = lsqnonlin( ...
    @(params) ODE_fit(params(1),params(2),params(3), ...
    t_stressor,stressor,t_em_intrp,emission_intrp), ...
    x0);

tsolv=[0 5];
ic = emission_intrp(1);
[t,sol] = ode45(@(t,g) ODE_eq(t,g,param_fit(1),param_fit(2), ...
    param_fit(3), t_stressor, stressor),tsolv,ic);
sol_intrp = interp1(t,sol,t_em_intrp);
plot(t_em_intrp,sol_intrp);