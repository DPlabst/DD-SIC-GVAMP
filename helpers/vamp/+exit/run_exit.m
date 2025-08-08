function [rateEXIT_end, nmseEXIT_end, rateEXIT_vec, nmseEXIT_vec, nuU1_smp, nuW1_smp] = run_exit( ...
        nu_vec, ...
        w_star, ...
        u_star, ...
        hnd_denoise_w1, ...
        hnd_denoise_u1, ...
        hnd_lmmse_u2w2, ...
        hnd_est_rates, ...
        idx_data, ...
        maxit, ...
        SNRlin, ...
        rGVAMP, ...
        rng_seed, ...
        EXIT_dbg_mode, ...
        Np, ...
        M ...
    )

    %% Create EXIT charts and find fixed point of GVAMP
    iterInf = inf; %Case where annealing in LMMSE denoiser is OFF
    lmmse_trig = 0; %Case where annealing in LMMSE denoiser is ON

    %% ---- Initialize ------
    nu_W1_vec = [];
    nu_W2_vec = [];
    nuW1_vs_nuW2 = [];
    nuW2_vs_nuW1 = []; 

    rng(rng_seed); %set seed

    % ------ Noise on U ------
    if rGVAMP %If real-valued prior
        Np1 = (1) * (randn(Np, 1)); %Fix realization
    else
        Np1 = (1 / sqrt(2)) * (randn(Np, 1) + 1j * randn(Np, 1)); %Fix realization
    end

    % ----- Noise on W ------
    Np2 = (1 / sqrt(2)) * (randn(M, 1) + 1j * randn(M, 1)); %Fix realization

    % Create handle for AWGN perturbed version of true symbols u_star
    u1 = @(nu_U1) u_star(idx_data) + sqrt(nu_U1) * Np1;

    if EXIT_dbg_mode
        se_curves_iter = [1, inf]; %Corresponds to annealing [ON, OFF]
    elseif ~EXIT_dbg_mode
        se_curves_iter = [inf]; %Corresponds to annealing [OFF]
    end

    for iter_val = se_curves_iter

        %% Handles that can also work with denoiser mismatch
        hnd_alpha1_vec = @(nu_U1_vec, nu_U1_mism_vec) exit.alpha1_vec(hnd_denoise_u1, u1, nu_U1_vec, nu_U1_mism_vec, rGVAMP); %Short handle that depends only on variances
        hnd_alphabeta2_vec = @(nu_U2_vec, nu_U2_mism_vec, nu_W2_vec, nu_W2_mism_vec) exit.alphabeta2_vec(hnd_lmmse_u2w2, u1, w_star, nu_U2_vec, nu_U2_mism_vec, nu_W2_vec, nu_W2_mism_vec, Np2, iter_val, lmmse_trig, rGVAMP); %Short handle that depends only on variances
        hnd_beta1_vec = @(nu_W1_vec, nu_W1_mism_vec) exit.beta1_vec(hnd_denoise_w1, w_star, nu_W1_vec, nu_W1_mism_vec, nu_vec, iter_val, Np2); %Short handle that depends only on variances

        if EXIT_dbg_mode %Generate transfer functions for EXIT chart

            nu_W1_vec = logspace(-5, 5, 500); % x-values
            nu_W2_vec = logspace(-5, 5, 500); % y-values

            [beta1_vs_nuW1] = hnd_beta1_vec(nu_W1_vec, []); % Output denoiser

            %% Transfer function: nuW2 = T1(nuW1)
            nuW2_vs_nuW1 = nu_W1_vec .* beta1_vs_nuW1 ./ (nu_W1_vec - beta1_vs_nuW1); % extrinsic

            %% Transfer function: nuW1 = T2(nuW2)
            nu_U2_vec = ones(1, length(nu_W1_vec)); %Dummy input, as LMMSE extrinsic (nuU1) is independent of nuU2

            [alpha2_vec, ~] = hnd_alphabeta2_vec(nu_U2_vec, [], nu_W2_vec, []); %LMMSE first part
            nu_U1_vec = nu_U2_vec .* alpha2_vec ./ (nu_U2_vec - alpha2_vec); %Convert to extrinsic

            alpha1_vec = hnd_alpha1_vec(nu_U1_vec, []); %Run input denoiser
            nuU2_vs_nuU1 = nu_U1_vec .* alpha1_vec ./ (nu_U1_vec - alpha1_vec); %Convert to extrinsic

            [~, beta2_vec] = hnd_alphabeta2_vec(nuU2_vs_nuU1, [], nu_W2_vec, []); %LMMSE second part
            nuW1_vs_nuW2 = nu_W2_vec .* beta2_vec ./ (nu_W2_vec - beta2_vec); %Convert to extrinsic
        end

        %% Calculate EXIT trajectory
        nuU2_smp = zeros(1, maxit);
        nuW2_smp = zeros(1, maxit);
        nuU1_smp = zeros(1, maxit);
        nuW1_smp = zeros(1, maxit);

        alpha1_smp = zeros(1, maxit);
        alpha2_smp = zeros(1, maxit);
        beta2_smp = zeros(1, maxit);

        rateEXIT_vec = zeros(1, maxit);
        nmseEXIT_vec = zeros(1, maxit);

        nuW1_smp(1) = 1e12; %Init
        nuW1_smp_mism(1) = nuW1_smp(1);
        nuW2_smp_mism = [];
        nuU1_smp_mism = [];
        nuU2_smp_mism = [];

        for k = 1:maxit
            if iter_val == inf %Matched EXIT trajectory

                %% First half-iteration Transfer function: nuW2 = T1(nuW1)
                beta1_ii = hnd_beta1_vec(nuW1_smp(k), []); %Matched output denoiser
                nuW2_smp(k) = nuW1_smp(k) * beta1_ii / (nuW1_smp(k) - beta1_ii); %EXT
                if nuW2_smp(k) <= 0; k = k - 1; break; end %Terminate

                %% Second half-iteration Transfer function: nuW1 = T2(nuW2)
                nu_U2_dum = ones(1); %Dummy input, as LMMSE extrinsic (nuU1) is independent of nuU2
                [alpha2_ii, ~] = hnd_alphabeta2_vec(nu_U2_dum, [], nuW2_smp(k), []); %Matched LMMSE first part
                alpha2_smp(k) = alpha2_ii + eps;

                if alpha2_smp(k) <= 0; k = k - 1; break; end %Terminate
                nuU1_smp(k) = nu_U2_dum .* alpha2_smp(k) ./ (nu_U2_dum - alpha2_smp(k)); %EXT

                [alpha1_ii, u1den] = hnd_alpha1_vec(nuU1_smp(k), []);
                alpha1_smp(k) = alpha1_ii + eps;
                if alpha1_smp(k) <= 0; k = k - 1; break; end %Terminate

                nuU2_smp(k) = nuU1_smp(k) .* alpha1_smp(k) / (nuU1_smp(k) - alpha1_smp(k)); %EXT

                [~, beta2_ii] = hnd_alphabeta2_vec(nuU2_smp(k), [], nuW2_smp(k), []); %LMMSE second part

                beta2_smp(k) = beta2_ii + eps;
                if beta2_smp(k) <= 0; k = k - 1; break; end %Terminate

                % Metrics:
                rateEXIT_vec(k) = hnd_est_rates(u1(nuU1_smp(k)).', nuU1_smp(k), rGVAMP);
                nmseEXIT_vec(k) = mean(abs(u_star(idx_data).' - u1den).^2) / SNRlin; %NMSE

                if k + 1 > maxit || (k > 1 && abs(rateEXIT_vec(k) - rateEXIT_vec(k - 1)) / abs(rateEXIT_vec(k)) < 1E-9)
                    break; %Terminate
                end

                nuW1_smp(k + 1) = nuW2_smp(k) .* beta2_smp(k) ./ (nuW2_smp(k) - beta2_smp(k)); %EXT

            else %Mismatched EXIT trajectory due to annealing

                %% First half-iteration Transfer function: nuW2 = T1(nuW1)
                [beta1, w1hat, p1] = hnd_beta1_vec(nuW1_smp(k), nuW1_smp_mism(k)); %Mismatched output denoiser
                beta1_smp(k) = mean(beta1);

                p2 = (nuW1_smp(k) * w1hat - beta1_smp(k) * p1.') / (nuW1_smp(k) - beta1_smp(k)); %What GVAMP believes
                nuW2_smp(k) = nuW1_smp(k) .* beta1_smp(k) ./ (nuW1_smp(k) - beta1_smp(k)); %EXT
                nuW2_smp_mism(k) = var(w_star.' - p2); %Actual MSE

                if nuW2_smp(k) <= 0; k = k - 1; break; end %Terminate

                %% Second half-iteration Transfer function: nuW1 = T2(nuW2)
                nu_U2_dum = ones(1); %Dummy input, as LMMSE extrinsic (nuU1) is independent of nuU2
                [alpha2_ii, ~, u2d_hat, ~, r2, ~] = hnd_alphabeta2_vec(nu_U2_dum, nu_U2_dum, nuW2_smp(k), nuW2_smp_mism(k)); %Matched LMMSE first part
                alpha2_smp(k) = alpha2_ii + eps;

                if alpha2_smp(k) <= 0; k = k - 1; break; end %Terminante
                r1 = (nu_U2_dum * u2d_hat - alpha2_smp(k) * r2) / (nu_U2_dum - alpha2_smp(k)); %EXT
                nuU1_smp(k) = nu_U2_dum .* alpha2_smp(k) ./ (nu_U2_dum - alpha2_smp(k)); %What GVAMP believes
                nuU1_smp_mism(k) = var(u_star(idx_data) - r1); %Actual MSE

                [alpha1_ii, u1den, r1] = hnd_alpha1_vec(nuU1_smp(k), nuU1_smp_mism(k)); %Run input denoiser
                alpha1_smp(k) = mean(alpha1_ii);

                if alpha1_smp(k) <= 0; k = k - 1; break; end

                r2 = (nuU1_smp(k) * u1den - alpha1_smp(k) * r1.') / (nuU1_smp(k) - alpha1_smp(k)); %EXT
                nuU2_smp(k) = nuU1_smp(k) .* alpha1_smp(k) / (nuU1_smp(k) - alpha1_smp(k)); %What GVAMP believes
                nuU2_smp_mism(k) = var(r2 - u_star(idx_data).'); %Actual MSE

                [~, beta2_ii, ~, w2_hat, ~, p2] = hnd_alphabeta2_vec(nuU2_smp(k), nuU2_smp_mism(k), nuW2_smp(k), nuW2_smp_mism(k)); %Matched LMMSE first part
                beta2_smp(k) = beta2_ii;

                if beta2_smp(k) <= 0; k = k - 1; break; end

                % Metrics:
                nmseEXIT_vec(k) = var(u_star(idx_data).' - u1den) / SNRlin; %Compute NMSE
                rateEXIT_vec(k) = hnd_est_rates(r1.', nuU1_smp_mism(k), rGVAMP);
                if k + 1 > maxit || (k > 1 && abs(rateEXIT_vec(k) - rateEXIT_vec(k - 1)) / abs(rateEXIT_vec(k)) < 1E-9)
                    break; %terminate
                end

                p1 = (nuW2_smp(k) * w2_hat - beta2_smp(k) * p2) / (nuW2_smp(k) - beta2_smp(k)); %EXT
                nuW1_smp(k + 1) = nuW2_smp(k) .* beta2_smp(k) ./ (nuW2_smp(k) - beta2_smp(k)); %What GVAMP believes
                nuW1_smp_mism(k + 1) = var(p1 - w_star); %Actual MSE
            end

        end

        rateEXIT_end = rateEXIT_vec(k);
        nmseEXIT_end = nmseEXIT_vec(k);

        % Some plotting
        plot_exit_dbg(EXIT_dbg_mode, iter_val, k, ...
            nuW1_smp, nuW2_smp, nuW1_smp_mism, nuW2_smp_mism, ...
            nu_W1_vec, nu_W2_vec, nuW1_vs_nuW2, nuW2_vs_nuW1, ...
            rateEXIT_vec, rateEXIT_end)

    end

end
