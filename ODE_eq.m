function dgdt = ODE_eq(t,g,w,c,k_d,st,s)
    s = interp1(st,s,t);
    dgdt = (max(g)./ (1+exp(-w*s + c)))-k_d*g;
end
