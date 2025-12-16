function [best_param,error_profile]=ODE_fit(starting_positions, p, t_plot,expected_emission,maximum)

    options = optimoptions(@lsqnonlin, ...
    'Algorithm','levenberg-marquardt');
    error_profile = zeros(1,length(starting_positions));
    best_param = zeros(length(starting_positions),3);
    for l=1:length(starting_positions)
        x0=repmat(starting_positions(l),1,3);

        [param_fit,resnorm] = lsqnonlin( ...
            @(params) ODE_solve(params(1),params(2),params(3), ...
                p,t_plot,expected_emission,maximum), ...
            x0, [0,0,0],[10,25,10],options);
        best_param(l,:) = param_fit;
        error_profile(l) = resnorm;
    end

end