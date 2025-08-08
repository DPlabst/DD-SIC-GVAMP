function [w1_hat, beta1_hat_vec, fk2, fk4, fk6] = denoise_w1(y, p1, nu_W1, nu, method)
    
    %% Output denoiser that computes messages (w1hat, beta1)
    % The output denoiser can consider models W = |W + N1|^2 + N2, where 
    % - N1 is _iid_ CSCG 
    % - N2 is _iid_ real Gaussian 
    % For most cases, N1 is treated in the LMMSE denoiser, so for the output denoiser var(N1) = 0

    nuN1 = nu(1); %(Mismatched) pre-intensity noise variance
    nuN2 = nu(2); %(Mismatched) post-intensity noise variance

    if nuN2 > 0 %If post-intensity noise
        nuWp = nu_W1 + nuN1;
        w = nuWp / 2;
        lambda = abs(p1).^2 / w;

        if strcmp(method, 'Davies') %May become numerically unstable for small variance 
            VPAflag = false;
            RelTol = 1E-13;
            AbsTol = 1E-13;

            fk2 = gx2_imhof_mod(y, w, 2, lambda, sqrt(nuN2), 0, RelTol, AbsTol, VPAflag); %Takes std of noise nuN2 as input
            fk4 = gx2_imhof_mod(y, w, 4, lambda, sqrt(nuN2), 0, RelTol, AbsTol, VPAflag); %Takes std of noise nuN2 as input
            fk6 = gx2_imhof_mod(y, w, 6, lambda, sqrt(nuN2), 0, RelTol, AbsTol, VPAflag); %Takes std of noise nuN2 as input

            fk4_o_fk2 = fk4 ./ (fk2 + eps);
            fk6_o_fk2 = fk6 ./ (fk2 + eps);

            abs2p = abs(p1).^2;

            if nuN1 == 0 %Simplifies, if nuN1 = 0
                w1_hat = p1 .* fk4_o_fk2;
                beta1_hat_vec = real((abs2p .* fk6_o_fk2 + nuWp * fk4_o_fk2) - abs(w1_hat).^2);

            else %General case
                w1_hat = (p1 / nuWp) .* (nuN1 + nu_W1 * fk4_o_fk2);
                beta1_hat_vec = real(nuN1 / nuWp^2 * (nuWp * nu_W1 + abs2p * nuN1) + ... %Part 1
                    (nu_W1^2 / nuWp^2) * (abs2p .* fk6_o_fk2 + nuWp * fk4_o_fk2) + ... %Part 2
                    2 * nuN1 * nu_W1 * (abs2p / nuWp^2) .* (fk4_o_fk2) + ... %Part 3
                    - abs(w1_hat).^2); %Part 4
            end

        elseif strcmp(method, 'VPA') % Very slow 

            syms a
            rho_arg = 2 * sqrt(a) .* abs(p1) ./ nuWp;
            f1 = besseli(0, rho_arg);
            f2 = besseli(2, rho_arg);
            abs2p = abs(p1).^2;

            pN2 = @(y) 1 / sqrt(2 * pi * nuN2) * exp(- (y - a).^2 ./ (2 * nuN2));
            E1 = exp(- (a + abs(p1).^2) / nuWp);

            fv2 = @(y) 1 / nuWp * f1 .* pN2(y) .* E1;
            fv4 = @(y) 1 ./ abs(p1) .* sqrt(a) * nuWp / (nuWp^2) .* besseli(1, rho_arg) .* pN2(y) .* E1;
            fv6 = @(y) a * 1 ./ (abs2p * nuWp) .* f2 .* pN2(y) .* E1;

            TolX = 1E-9;
            fk2 = double(vpaintegral(fv2(y), 0, inf, 'RelTol', TolX, 'AbsTol', TolX)); %Correct
            fk4 = double(vpaintegral(fv4(y), 0, inf, 'RelTol', TolX, 'AbsTol', TolX));
            fk6 = double(vpaintegral(fv6(y), 0, inf, 'RelTol', TolX, 'AbsTol', TolX));

            fk4_o_fk2 = fk4 ./ (fk2);
            fk6_o_fk2 = fk6 ./ (fk2);

            w1_hat = (p1 / nuWp) .* (nuN1 + nu_W1 * fk4_o_fk2);

            beta1_hat_vec = real(nuN1 / nuWp^2 * (nuWp * nu_W1 + abs2p * nuN1) + ... %Part 1
                (nu_W1^2 / nuWp^2) * (abs2p .* fk6_o_fk2 + (nuWp * fk4_o_fk2)) + ... %Part 2
                2 * nuN1 * nu_W1 * (abs2p / nuWp^2) .* (fk4_o_fk2) + ... %Part 3
                - abs(w1_hat).^2); %Part 4


        elseif strcmp(method, 'SPA')

            [fk2, fk2_pre, fk2_logE] = wdenoise.spa_pdf_vec(y, nuWp, p1, 2, nuN2); %variance of noise nuN2 as input
            [fk4, fk4_pre, fk4_logE] = wdenoise.spa_pdf_vec(y, nuWp, p1, 4, nuN2); %variance of noise nuN2 as input
            [fk6, fk6_pre, fk6_logE] = wdenoise.spa_pdf_vec(y, nuWp, p1, 6, nuN2); %variance of noise nuN2 as input

            fk4_o_fk2 = fk4_pre ./ fk2_pre .* exp(fk4_logE - fk2_logE);
            fk6_o_fk2 = fk6_pre ./ fk2_pre .* exp(fk6_logE - fk2_logE);

            abs2p = abs(p1).^2;

            if nuN1 == 0 % Simplifies
                w1_hat = p1 .* fk4_o_fk2;
                beta1_hat_vec = real((abs2p .* fk6_o_fk2 + nuWp * fk4_o_fk2)) - abs(w1_hat).^2;

            else %General case
                w1_hat = (p1 / nuWp) .* (nuN1 + nu_W1 * fk4_o_fk2);
                beta1_hat_vec = real(nuN1 / nuWp^2 * (nuWp * nu_W1 + abs2p * nuN1) + ... %Part 1
                    (nu_W1^2 / nuWp^2) * (abs2p .* fk6_o_fk2 + nuWp * fk4_o_fk2) + ... %Part 2
                    2 * nuN1 * nu_W1 * (abs2p / nuWp^2) .* (fk4_o_fk2) + ... %Part 3
                    - abs(w1_hat).^2); %Part 4
            end
        end

    elseif nuN1 >= 0 && nuN2 == 0 %If pre-intensity noise, closed-form solutions are available

        [I1overI0, rho_arg] = wdenoise.I1_div_I0(y, p1, nu_W1, nuN1);

        w1_hat = p1 ./ (nu_W1 / nuN1 + 1) + ...
            sqrt(y) .* exp(1j * angle(p1)) ./ (nuN1 ./ nu_W1 + 1) .* ...
            I1overI0;

        beta1_hat_vec = y ./ (1 + nuN1 ./ nu_W1).^2 + abs(p1).^2 ./ (1 + nu_W1 / nuN1).^2 + ...
            (1 + rho_arg .* I1overI0) ./ (1 / nuN1 + 1 ./ nu_W1) - abs(w1_hat).^2; %Subtract abs(CME)^2

    end

end
