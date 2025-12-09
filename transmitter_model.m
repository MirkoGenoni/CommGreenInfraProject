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


t_new= linspace(0,30,10000*39);
y = param_fit(1)+param_fit(2).*t_new+param_fit(3).*t_new.^2+ ...
    param_fit(4).*t_new.^3 + param_fit(5).*t_new.^4;

figure;
scatter(data_stressor.x, data_stressor.EmissionToLeafArea);
hold on; 
plot(t_new, y);
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
sampling_T = 1;
sampling_f = 1/sampling_T;
interpolated_length = (sampling_f*length(data_emission.x))-1;
new_x = zeros(1, interpolated_length);
interpolated_data = zeros(1,interpolated_length);

for j=2:length(data_emission.x)
    curr_interval = (j-2)*(sampling_f)+1:(j-1)*(sampling_f)+1;

    new_x(1,curr_interval) = ...
        linspace(data_emission.x(j-1),data_emission.x(j),(sampling_f)+1);

    interpolated_data(1,curr_interval)= ...
        linspace(data_emission.larvae_2(j-1),data_emission.larvae_2(j), ...
        (sampling_f)+1);
end

% EXTRACT STRESSOR FROM CURRENT EMISSION
larvae_2 = data_emission.larvae_2';
leaf_cons = zeros(1,interpolated_length);
for index1= 1:length(interpolated_data)
    temp = find(round(y,4)==round(interpolated_data(index1),4));
    leaf_cons(1,index1) = y(round((max(temp) + min(temp))/2)-1)*1e-4;
end

% LOSS FUNCTION BETWEEN OUR EQUATION AND THE EXPERIMENTAL DATA
delta_value = interpolated_data(2:end)-interpolated_data(1:end-1);
equation = @(parameters) ((max(interpolated_data)./ ... 
    (1+exp(-parameters(1)*leaf_cons(1:end-1) + parameters(2)))) ...
    -parameters(3)*interpolated_data(1:end-1))- ...
    delta_value;

% FITTING
x0=[1e-1,1e-1,1e-1]; %initial parameters
options = optimoptions('lsqnonlin','Display','iter');
[param_fit,resnorm,residual,eflag,output] = lsqnonlin(equation,x0, ...
    [], [], options);

% CONSTRUCT FUNCTION FROM PARAMETERS FROM THE FITTING
fitted = (max(interpolated_data)./ ... 
    (1+exp(-param_fit(1)*leaf_cons + param_fit(2)))) ...
    -param_fit(3)*interpolated_data;


estimated = zeros(1, interpolated_length);
for index = 1:length(interpolated_data)
    if index == 1
        estimated(index) = interpolated_data(1);
    else
        estimated(index) = estimated(index-1)+fitted(index-1);
    end
end

original_var = sum((data_emission.larvae_2 - ...
    mean(data_emission.larvae_2)).^2);
error = sum((estimated-interpolated_data).^2);
r_2_esmission = 1-(error/original_var)

figure;
plot(new_x,interpolated_data);
hold on;
scatter(data_emission.x, data_emission.larvae_2);
plot(new_x,estimated);
hold off;
s.Marker="square";
s.MarkerFaceColor="red";
fontsize(16,"points");
xlabel("Days"); ylabel("LOX emission rate [nmol m^{-2} s^{-1}]", ...
    'Interpreter', 'tex');
grid on;