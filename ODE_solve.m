function out = ODE_solve(w,c,k_d, p, t_plot,expected_emission,maximum)
    tsolv=[0 5];
    ic = expected_emission(1);
    sol = ode45(@(t,g) ODE_eq(t,g,w,c,k_d,p,maximum), ...
        tsolv,ic);

    solution_intrp = deval(sol,t_plot);
    out = solution_intrp - expected_emission;
end