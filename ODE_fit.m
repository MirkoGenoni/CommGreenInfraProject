function [best_param,error_profile]=ODE_fit(starting_positions, p, t_plot,expected_emission)

    options = optimoptions(@lsqnonlin, ...
    'Algorithm','levenberg-marquardt');
    options.MaxIterations = 1000;
    error_profile = zeros(1,length(starting_positions));
    prev_resnorm = 1000;
    for l=1:length(starting_positions)
        x0=repmat(starting_positions(l),1,3);

        [param_fit,resnorm,residual,eflag,output] = lsqnonlin( ...
            @(params) ODE_solve(params(1),params(2),params(3), ...
                p,t_plot,expected_emission), ...
            x0, [],[],options);

        error_profile(l) = resnorm;
        if(prev_resnorm>resnorm)
            best_param=param_fit;
        end
        prev_resnorm = resnorm;
    end

end