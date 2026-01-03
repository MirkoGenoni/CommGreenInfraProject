function C_air = anisotropic_gaussian_puff(Q,class,t,x_rx,y_rx,z_rx,h,M,delta_p,eps)

    % Parameters setup
    u = class.u;
    r_y = sqrt(class.sigma_y)/2;
    r_z = sqrt(class.sigma_z)/2;
    C_air = 0;

    for i = 1:length(y_rx) %h more than 1 element it's a MISO channel
        % Concentration in air
        for j = 1:1:M;    
            t = j*delta_p;
            A = Q(j)./(8.*(pi.*u.*sqrt(r_y.*r_z)).^(3/2));
            B = exp(-((x_rx-u.*t).^2)./(4.*r_y));
            C = exp(-(y_rx(i).^2)./(4.*r_y));
            D = exp(-((z_rx-h).^2)./(4.*r_z))+exp(-((z_rx+h).^2)./(4.*r_z));
            C_new = A.*B.*C.*D;

            if (C_new<eps)
                C_new = 0;
            end
            
            C_air = C_air + C_new;

    end

end