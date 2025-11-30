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
    data_stressor.Curve1);

a = param_fit(1);
b = param_fit(2);

x= 0:0.01:30;
y=a*exp(b*x);

figure;
scatter(data_stressor.x, data_stressor.Curve1);
hold on;
plot(x, y);

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
larvae_2 = data_emission.larvae_2;

leaf_cons = zeros(length(larvae_2),1);
for index1= 1:length(larvae_2)
    temp = find(round(y,4)==round(larvae_2(index1),4));
    leaf_cons(index1,1) = (round((max(temp) + min(temp))/2)-1)*0.01;
end

equation = @(parameters, z) (max(larvae_2)./ ... 
    (1+exp(-parameters(1)*leaf_cons(z) + parameters(2)))) ...
    -parameters(3)*larvae_2(z);

I_measured = larvae_2(2:end)-larvae_2(1:end-1);
I_measured = [I_measured;0];

initial_guess = [10,24,8.5];

options = optimoptions('lsqcurvefit', ...
    'Algorithm','levenberg-marquardt', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-12, ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 10000, ...
    'Display','iter-detailed');

[param_fit, r] = lsqcurvefit(equation, initial_guess, x+1, ...
    I_measured, [], [], options);

w_fit = param_fit(1);
c_fit = param_fit(2);
k_d_fit = param_fit(3);

I_est = (max(larvae_2)/ ... 
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

%--------------------------------------------------------------------------