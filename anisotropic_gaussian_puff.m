function C_air = anisotropic_gaussian_puff(Q,class,t,x_rx,y_rx,z_rx,h)

    % Parameters setup
    u = class.u;
    r_y = sqrt(class.sigma_y)/2;
    r_z = sqrt(class.sigma_z)/2;
    C_air = 0;

    for i = 1:length(y_rx) %h more than 1 element it's a MISO channel
        % Concentration in air
        A = Q./(8.*(pi.*u.*sqrt(r_y.*r_z)).^(3/2));
        B = exp(-((x_rx-u.*t).^2)./(4.*r_y));
        C = exp(-(y_rx(i).^2)./(4.*r_y));
        D = exp(-((z_rx-h).^2)./(4.*r_z))+exp(-((z_rx+h).^2)./(4.*r_z));
    
        C_air = C_air + A.*B.*C.*D;

    end

end