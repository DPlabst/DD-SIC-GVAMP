function [logpYgZ] = log_pygw(y, w, nu)

    % Compute p(y|w)
    nuN1 = nu(1); 
    nuN2 = nu(2); 

    if nuN1 == 0
        nuN1 = nuN1 + 1E-6; 
    end 
    if nuN2 > 0 %If post-intensity noise 
        if nuN1 > 0 %If also pre-intensity noise 
            nuWp = nuN1;
            [~, fk2_pre, fk2_logE] = wdenoise.spa_pdf_vec(y, nuWp, w, 2, nuN2);
            logpYgZ = log(fk2_pre) + fk2_logE;
        else
            logpYgZ = log(1 ./ sqrt(2 * pi * nuN2)) - abs(y - abs(w).^2).^2 / (2 * nuN2); %Conditionally Gaussian
        end

    elseif nuN1 > 0 && nuN2 == 0 %If only pre-intensity noise: closed-form  available

        rho_arg_lkh = 2 * sqrt(y) .* abs(w) ./ nuN1; 
        bessel_f_sc_I0_lkh = besseli(0, rho_arg_lkh, 1);

        logpYgZ = log(double(y >= 0)) + log(1 / nuN1) ...
            - (y + abs(w).^2) / nuN1 + log(bessel_f_sc_I0_lkh) ...
            + (abs(real(rho_arg_lkh)));

    end
