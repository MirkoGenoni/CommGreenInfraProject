function dgdt = ODE_eq(t,g,w,c,k_d,p,maximum)
    s = polyval(p,t);
    dgdt = maximum./(1+exp(-w*s + c))-k_d*g;
end
