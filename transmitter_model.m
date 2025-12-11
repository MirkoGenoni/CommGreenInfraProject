%%
close all; clc; clear;

<<<<<<< HEAD
%% ----- STRESSOR FITTING --------------------------------------------------
=======
% ---------- LEAF AREA FITTING --------------------------------------------
>>>>>>> b92f287f5abd4dc944a6995ae5f32fde5c9d511c

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
<<<<<<< HEAD
scatter(data_stressor.x, data_stressor.Curve1);
hold on;
plot(x, y);
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");
fontsize(16,"points");

% -------------------------------------------------------------------------

%% ----- EMISSION FITTING --------------------------------------------------

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
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10, ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 100000, ...
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
=======
scatter(data_stressor.x, data_stressor.EmissionToLeafArea);
hold on; 
plot(t_stressor, stressor);
>>>>>>> b92f287f5abd4dc944a6995ae5f32fde5c9d511c
grid on;
xlabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
ylabel("LOX emission [nmol m^{-2} s^{-1}]");
fontsize(16,"points");

%--------------------------------------------------------------------------


% --------- EMISSION FITTING ----------------------------------------------

<<<<<<< HEAD
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

%%I_actual = data_emission_LS.Curve1(2:end)-data_emission_LS.Curve1(1:end-1);
I_actual = larvae_2(2:end)-larvae_2(1:end-1);

max = 0.05;
w = 10.02;
c = 23.83;
k_d = 8.52;

curr_s = (c-log(((max./(I_actual+k_d.*data_emission_LS.Curve1(1:end-1))))-1))./w;

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
=======
% DATA SANIFICATION
>>>>>>> b92f287f5abd4dc944a6995ae5f32fde5c9d511c
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

%PLOT STRESSOR PROFILE
figure;
plot(t_em_intrp, leaf_cons);
xlabel("Days");
ylabel("Leaf area eaten [cm^2]", 'Interpreter', 'tex');
fontsize(16,"points");

figure;
plot(t_em_intrp,emission_intrp);
hold on;
scatter(t,emission);

% LEAST SQUARE ON DIFFERENTIAL EQUATION
x0=[1,1,1];
[param_fit,resnorm,residual,eflag,output] = lsqnonlin( ...
    @(params) ODE_fit(params(1),params(2),params(3), ...
    t_em_intrp,leaf_cons,t_em_intrp,emission_intrp), ...
    x0);

% SOLVE DIFFERENTIAL EQUATION WITH OPTIMAL PARAMETERS
tsolv=[0 5];
ic = emission_intrp(1);
[t,sol] = ode45(@(t,g) ODE_eq(t,g,param_fit(1),param_fit(2), ...
    param_fit(3), t_em_intrp, leaf_cons),tsolv,ic);
sol_intrp = interp1(t,sol,t_em_intrp);
plot(t_em_intrp,sol_intrp);
xlabel("Days");
ylabel("LOX emission [nmol m^{-2} s^{-1}]", 'Interpreter', 'tex');
fontsize(16,"points");
legend('Emission interpolation', 'Experimental points', ...
    'Emission fitted to differential eq.');
