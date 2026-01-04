%%
close all; close all hidden; clc; clear;

% --------- EMISSION FITTING ----------------------------------------------
larve_num = 4;
show = "best";
show_results = 6;
dt= 1;
polynomial_grade = 1:2;
starting_positions= [linspace(10^-3,10^-2,10), ...
    linspace(2*10^-2,10^-1,9), linspace(2*10^-1,1,9), ...
    linspace(2,10,9)]; %STARTING VALUES FOR THE PARAMETERS

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

emission_profiles = [data_emission.larvae_2,data_emission.larvae_4, ...
    data_emission.larvae_8];

% EMISSION INTERPOLATION
t=0:5;
index=log2(larve_num);
emission=emission_profiles(:,index);
maximum=max(emission_profiles(:, index));

t_em_intrp=0:dt:5;

emission_intrp = interp1(t, emission,t_em_intrp,'linear');
original_var_emission = sum((emission_intrp-mean(emission_intrp)).^2);

% STRESSOR AS NUMBER OF LARVAE
leaf_cons = zeros(1,length(t_em_intrp));
leaf_cons(1:end-1) = larve_num;
leaf_cons(end) = 0;

for k=1:length(polynomial_grade)
    % FITTING POLYNOMIAL TO STRESS PROFILE
    r_2_leaf_area = zeros(1, length(polynomial_grade));
    original_var_stress = sum((leaf_cons - mean(leaf_cons)).^2);
    p = polyfit(t_em_intrp, leaf_cons, polynomial_grade(k));
    leaf_fit = polyval(p, t_em_intrp);
    residuals = leaf_cons - leaf_fit;
    resnorm_sq = sum(residuals.^2);
    r_2_leaf_area(k) = 1-(resnorm_sq/original_var_stress);

    %PLOT STRESSOR PROFILE
    figure;
    scatter(t_em_intrp, leaf_cons);
    hold on;
    plot(t_em_intrp, leaf_fit);
    xlabel("Days");
    ylabel("Number of larvae");
    fontsize(16,"points");
    
    %PLOT EXPERIMENTAL EMISSION POINTS
    figure;
    plot(t_em_intrp,emission_intrp);
    hold on;
    scatter(t,emission);
    
    % LEAST SQUARE ON DIFFERENTIAL EQUATION
    [fit_param, error_profile] = ODE_fit(starting_positions,p, ...
        t_em_intrp,emission_intrp,maximum, original_var_emission);
    
    %EXTRACT ONLY PARAMETERS WITH THE BEST r^2
    if show == "best"
        maximum_r_2 = maxk(error_profile,show_results);
        indexes_best_values = zeros(show_results,1);
        for j=1:length(maximum_r_2)
            indexes_best_values(j) = find(error_profile==maximum_r_2(j));
        end
        indexes_best_values = sort(indexes_best_values,1,"ascend");
        fit_param = fit_param(indexes_best_values,:);
        error_profile = error_profile(indexes_best_values);
        starting_positions=starting_positions(indexes_best_values);
    end

    %PRINT TABLE WITH PARAMETERS, STARTING CONDITIONS AND r^2
    fit_param_table = [starting_positions', fit_param, error_profile'];
    fit_param_table = array2table(fit_param_table);
    fit_param_table.Properties.VariableNames(1:5) = {'Initial postion', ...
        'w', 'c', 'k_d', 'r_2'};
    fig = uifigure(Name=strcat("POLY DEGREE: ", ...
        num2str(polynomial_grade(k))));
    uit = uitable(fig,"Data",fit_param_table);

    % SOLVE DIFFERENTIAL EQUATION WITH OPTIMAL PARAMETERS
    tsolv=[0 5];
    ic = emission_intrp(1);
    for index=1:size(fit_param,1)
        [t_solver,sol] = ode45(@(t,g) ODE_eq(t,g,fit_param(index,1), ...
            fit_param(index,2), fit_param(index,3), p,maximum),tsolv,ic);
        % PLOT FITTED SOLUTION TO DIFFERENTIAL EQUATION
        plot(t_solver,sol);
        xlabel("Days");
        ylabel("LOX emission [nmol m^{-2} s^{-1}]", 'Interpreter', 'tex');
        fontsize(16,"points");
    end
    labels = [{'Emission interpolation', 'Experimental points'}, ...
              compose('%g', starting_positions)];
    
    legend(labels)
    
    %PLOT ERROR PROFILE OF FITTED DIFF. EQ.
    figure;
    scatter(starting_positions,error_profile);
    xlabel("Starting positions");
    ylabel('r^2', 'Interpreter','tex');
    fontsize(16,"points");
end
%%
function dgdt = ODE_eq(t,g,w,c,k_d,p,maximum)
    s = polyval(p,t);
    dgdt = maximum./(1+exp(-w*s + c))-k_d*g;
end

function out = ODE_solve(w,c,k_d, p, t_plot,expected_emission,maximum)
    tsolv=[0 5];
    ic = expected_emission(1);
    sol = ode45(@(t,g) ODE_eq(t,g,w,c,k_d,p,maximum), ...
        tsolv,ic);

    solution_intrp = deval(sol,t_plot);
    out = solution_intrp - expected_emission;
end

function [fit_param,error_profile]=ODE_fit(starting_positions, p, ...
    t_plot,expected_emission,maximum, original_var)

    options = optimoptions(@lsqnonlin, ...
    'Algorithm','levenberg-marquardt');
    error_profile = zeros(1,length(starting_positions));
    fit_param = zeros(length(starting_positions),3);
    for l=1:length(starting_positions)
        x0=repmat(starting_positions(l),1,3);

        [param_fit,resnorm] = lsqnonlin( ...
            @(params) ODE_solve(params(1),params(2),params(3), ...
                p,t_plot,expected_emission,maximum), ...
            x0, [0,-50,-50],[50,50,50],options);
        fit_param(l,:) = param_fit;
        error_profile(l) = 1-(resnorm/original_var);
    end

end