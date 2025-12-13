function dgdt = ODE_eq(t,g,w,c,k_d,p)
    s = polyval(p,t);
    dgdt = max(g)./(1+exp(-w*s + c))-k_d*g;
end
