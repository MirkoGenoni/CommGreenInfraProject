clc; clear; close all;

% ----- STRESSOR FITTING --------------------------------------------------

file = fopen("Data/Transmitter/LOX_to_leaf.csv", "r");
sanizited_data = ...
    fopen("Data/Transmitter/sanitized_LOX_to_leaf_values.csv", "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
    disp(new_line);
end
fclose(file);
fclose(sanizited_data);


data_stressor = ...
    readtable("Data/Transmitter/sanitized_LOX_to_leaf_values.csv", ...
    'PreserveVariableNames',true);

graph = @(param, x) param(1)*exp(param(2)*x);

initial_guess = [1e-1,1e-1];
param_fit = lsqcurvefit(graph, initial_guess, data_stressor.x, ...
    data_stressor.EmissionToLeafArea);

a = param_fit(1);
b = param_fit(2);

x= 0:0.01:30;
y=a*exp(b*x);

figure;
scatter(data_stressor.x, data_stressor.EmissionToLeafArea);
hold on;
plot(x, y);
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");
fontsize(16,"points");

% -------------------------------------------------------------------------

% ----- EMISSION FITTING --------------------------------------------------

file = fopen("Data/Transmitter/emission.csv", "r");
sanizited_data = fopen("Data/Transmitter/sanitized_emission_data.csv", ...
    "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
    disp(new_line);
end
fclose(file);
fclose(sanizited_data);

data_emission = readtable("Data/Transmitter/sanitized_emission_data.csv", ...
    'PreserveVariableNames',true);

x = 0:5;
larvae_2 = data_emission.larvae_2';

leaf_cons = zeros(1,length(larvae_2));
for index1= 1:length(larvae_2)
    temp = find(round(y,4)==round(larvae_2(index1),4));
    leaf_cons(1,index1) = (round((max(temp) + min(temp))/2)-1)*0.01*1e-4;
end

equation = @(parameters, z) (max(larvae_2)./ ... 
    (1+exp(-parameters(1)*leaf_cons(z+1) + parameters(2)))) ...
    -parameters(3)*larvae_2(z+1);

I_measured = larvae_2(2:end)-larvae_2(1:end-1);

initial_guess = [0.1,0.1,0.1];

options = optimoptions('lsqcurvefit', ...
    'Algorithm','levenberg-marquardt', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 10000, ...
    'Display','iter-detailed');

[param_fit, r] = lsqcurvefit(equation, initial_guess, x(1:end-1), ...
    I_measured(1:end), [], [], options);

w_fit = param_fit(1);
c_fit = param_fit(2);
k_d_fit = param_fit(3);

I_est = (max(larvae_2)./ ... 
    (1+exp(-w_fit*leaf_cons(1:end-1) + c_fit))) ...
    -k_d_fit.*larvae_2(1:end-1);


estimated = zeros(1, length(I_est));
for index = 1:length(x)
    if index == 1
        estimated(index) = larvae_2(1);
    else
        estimated(index) = estimated(index-1)+I_est(index-1);
    end
end

figure;
scatter(x, larvae_2);
hold on;
plot(x, estimated);
grid on;
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");

fontsize(16,"points");

%--------------------------------------------------------------------------

%% OBTAIN STRESSOR FROM THEIR DATA

file = fopen("Data/Transmitter/LOX_emission_LS.csv", "r");
sanizited_data = fopen("Data/Transmitter/sanitized_LOX_emission_LS.csv", ...
    "w");

while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
    disp(new_line);
end
fclose(file);
fclose(sanizited_data);

data_emission_LS = readtable("Data/Transmitter/sanitized_LOX_emission_LS.csv", ...
    'PreserveVariableNames',true);

I_actual = data_emission_LS.EmissionToLeafArea(2:end)-data_emission_LS.EmissionToLeafArea(1:end-1);

max = 0.05;
w = 10.02;
c = 23.83;
k_d = 8.52;

curr_s = (c-log(((max./(I_actual+k_d.*data_emission_LS.EmissionToLeafArea(1:end-1))))-1))./w;

% OBTAIN MEASURE FROM THEIR PARAMETERS

s = 2.5;
max = 0.05;
w = 10.02;
c = 23.83;
k_d = 8.52;

I_est = (max./ ... 
    (1+exp(-w*s + c))) ...
    -k_d.*larvae_2(1:end-1);
estimated = zeros(1, length(I_est));

for index = 1:length(x)
    if index == 1
        estimated(index) = larvae_2(1);
    else
        estimated(index) = estimated(index-1)+I_est(index-1);
    end
end

figure;
scatter(x, larvae_2);
hold on;
plot(x, estimated);
grid on;
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");
fontsize(16,"points");

%% PLOT POINTS
close all; clc; clear;
file = fopen("Data/Transmitter/emission.csv", "r");
sanizited_data = fopen("Data/Transmitter/sanitized_emission_data.csv", ...
    "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
    disp(new_line);
end
fclose(file);
fclose(sanizited_data);

data_emission = readtable("Data/Transmitter/sanitized_emission_data.csv", ...
    'PreserveVariableNames',true);


figure;
s= scatter(data_emission.x, data_emission.larvae_2);
s.Marker="square";
s.MarkerFaceColor="red";
fontsize(16,"points");
xlabel("Days"); ylabel("LOX emission rate [nmol m^{-2} s^{-1}]", ...
    'Interpreter', 'tex');
grid on;


figure;
s=scatter(data_emission.x, data_emission.larvae_4);
s.Marker="square";
s.MarkerFaceColor="red";
fontsize(16,"points");
xlabel("Days"); ylabel("LOX emission rate [nmol m^{-2} s^{-1}]", ...
    'Interpreter', 'tex');
grid on;


figure;
s=scatter(data_emission.x, data_emission.larvae_8);
s.Marker="square";
s.MarkerFaceColor="red";
fontsize(16,"points");
xlabel("Days"); ylabel("LOX emission rate [nmol m^{-2} s^{-1}]", ...
    'Interpreter', 'tex');
grid on;

%% PLOT AREA LEAF EATEN
close all; clc; clear;

file = fopen("Data/Transmitter/LOX_to_leaf.csv", "r");
sanizited_data = ...
    fopen("Data/Transmitter/sanitized_LOX_to_leaf_values.csv", "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
    disp(new_line);
end
fclose(file);
fclose(sanizited_data);


data_stressor = ...
    readtable("Data/Transmitter/sanitized_LOX_to_leaf_values.csv", ...
    'PreserveVariableNames',true);

graph = @(param, x) param(1)+param(2).*x+param(3).*x.^2+param(4).*x.^3+ ... 
    param(5).*x.^4;

initial_guess = [1e-1,1e-1,1e-1,1e-1,1e-1,1e-1];

param_fit = lsqcurvefit(graph, initial_guess, data_stressor.x, ...
    data_stressor.EmissionToLeafArea);

x= 0:0.01:30;
y=param_fit(1)+param_fit(2).*x + param_fit(3).* x.^2+param_fit(4).*x.^3+ ...
param_fit(5).*x.^4;

figure;
scatter(data_stressor.x, data_stressor.EmissionToLeafArea);
hold on;
plot(x, y);
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");

%%
close all; clc; clear;

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

original_var = sum((data_stressor.EmissionToLeafArea-mean(data_stressor.EmissionToLeafArea)).^2);

r_2 = 1-(resnorm/original_var);


t_new= linspace(0,30,10000*39);
y = param_fit(1)+param_fit(2).*t_new+param_fit(3).*t_new.^2+param_fit(4).*t_new.^3 ...
    + param_fit(5).*t_new.^4;
figure;
scatter(data_stressor.x, data_stressor.EmissionToLeafArea);
hold on; 
plot(t_new, y);

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
interpolated_data = zeros(1,(1/sampling_T)*(length(data_emission.x)-1));
for j=2:length(data_emission.x)
    new_x(1,(j-2)*(1/sampling_T)+1:(j-1)*(1/sampling_T)+1) = ...
        linspace(data_emission.x(j-1),data_emission.x(j),(1/sampling_T)+1);
    interpolated_data(1,(j-2)*(1/sampling_T)+1:(j-1)*(1/sampling_T)+1)= ...
        linspace(data_emission.larvae_2(j-1),data_emission.larvae_2(j), ...
        (1/sampling_T)+1);
end

% FITTING
larvae_2 = data_emission.larvae_2';
leaf_cons = zeros(1,length(interpolated_data));
for index1= 1:length(interpolated_data)
    temp = find(round(y,4)==round(interpolated_data(index1),4));
    leaf_cons(1,index1) = (round((max(temp) + min(temp))/2)-1)*0.01*1e-4;
end

equation = @(parameters) ((max(interpolated_data)./ ... 
    (1+exp(-parameters(1)*leaf_cons(1:end-1) + parameters(2)))) ...
    -parameters(3)*interpolated_data(1:end-1))- ...
    (interpolated_data(2:end)-interpolated_data(1:end-1));


x0=[1e-1,1e-1,1e-1];

options = optimoptions('lsqnonlin','Display','iter');
[param_fit,resnorm,residual,eflag,output] = lsqnonlin(equation,x0, [], [], options);

fitted = (max(interpolated_data)./ ... 
    (1+exp(-param_fit(1)*leaf_cons + param_fit(2)))) ...
    -param_fit(3)*interpolated_data;

for index = 1:length(interpolated_data)
    if index == 1
        estimated(index) = interpolated_data(1);
    else
        estimated(index) = estimated(index-1)+fitted(index-1);
    end
end

original_var = sum((data_emission.larvae_2-mean(data_emission.larvae_2)).^2);
error = sum((estimated-interpolated_data).^2);
r_2 = 1-(error/original_var);

figure;
plot(new_x,interpolated_data);
hold on;
scatter(data_emission.x, data_emission.larvae_2);
plot(new_x,estimated);
hold off;
