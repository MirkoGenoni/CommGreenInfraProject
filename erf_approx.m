function y = erf_approx(x, N)

    y = 0;

    for n = 0:N
        term = (-1)^n.*x.^(2*n+1)./(factorial(n)*(2*n+1));
        y = y+term;
    end

    y = (2/sqrt(pi))*y
end