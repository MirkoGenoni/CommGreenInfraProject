function [C_air, C_air_far_field] = anisotropic_gaussian_puff(Q,class,t,x_rx,y_rx,z_rx,h,M,delta_p,eps)

    % Parameters setup
    u = class.u;
    r_y = (class.sigma_y).^2/2;
    r_z = (class.sigma_z).^2/2;

    C_air = 0;
    C_air_far_field = 0;

    count = 0;
    max_count = length(y_rx);

    for i = 1:length(y_rx) %h more than 1 element it's a MISO channel
        % Concentration in air
        
        parfor j = 1:1:round(M) %togli il far field
            t = j*delta_p;

            % A = Q(round(M)+1-j)./(8.*(pi.*u.*sqrt(r_y.*r_z)).^(3/2));
            % B = exp(-((x_rx-u.*t).^2)./(4.*r_y));
            % C = exp(-(y_rx(i).^2)./(4.*r_y));
            % D = exp(-((z_rx-h).^2)./(4.*r_z))+exp(-((z_rx+h).^2)./(4.*r_z));
            % C_new = A.*B.*C.*D;
            
            A_sigma = Q(round(M)+1-j)./((class.sigma_y.*class.sigma_y.*class.sigma_z).*(2*pi).^(3/2));
            B_sigma = exp(-((x_rx-u.*t).^2)./(2*class.sigma_y.^2));
            C_sigma = exp(-(y_rx(i).^2)./(2*class.sigma_y.^2));
            D_sigma = exp(-((z_rx-h).^2)./(2*class.sigma_z.^2))+exp(-((z_rx+h).^2)./(2*class.sigma_z.^2));
            C_air_sigma = A_sigma.*B_sigma.*C_sigma.*D_sigma;

            % if (C_new<eps)
            %     C_new = 0;
            % end

            if (C_air_sigma<eps)
                C_air_sigma = 0;
            end
            
            C_air = C_air + C_air_sigma;
            % C_air_far_field = C_air_far_field + C_new;

        end

            %count = count+1;
            %display((count/max_count)*100,'Completed');

    end

end