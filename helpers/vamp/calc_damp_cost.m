function [J] = calc_damp_cost(y, ...
        N, ...
        r1, nu_U1, ... %VAMP belief mean and variance
        u1, alpha1, ... %VAMP belief mean and variance
        nu_N1sim, ... %Pre-intensity noise on whole simulation bandwidth [-B,+B]
        gamma, ... %Bd/Bsam
        nu_p_vec, ... % Assumed noise [Pre, Post]
        pU, ... %Input PMF
        Adata, ...
        z_intf, ...
        trAdAdh, ...
        idx_data, ...
        Hcfft, ...
        Hpfft, ...
        Pcfg, ...
        N_os,...
        rGVAMP, ...
        hnd_denoise_u1 ... 
    )

    %% Cost function for adaptive damping 
    M = N_os * N; %Number of observations 
    [~, ~, log_Q_UgR] = hnd_denoise_u1(r1, nu_U1, rGVAMP); %Approximate posterior belief 
    DivQP = log2(exp(1)) * mean(sum(exp(log_Q_UgR) .* (log_Q_UgR - log(pU.')), 1));

    %Approximate cost EQ_log2pYgU := int_u Q(u|y) * log2(p(y|v^{\ell-1},u))
    if any(contains(["cDC", "cDR", "cU", "cO", "_"], Pcfg)) %FFT-based 
        u2d_hat_padded = zeros(N, 1);
        u2d_hat_padded(idx_data) = u1; %T*x
        u2d_t1 = Hpfft .* (1/sqrt(N)*fft(u2d_hat_padded, N)); %Note: unitary n-fft 
        u2d_t2 = sqrt(2*N)/sqrt(N_os) * ifft(Hcfft .* [u2d_t1; u2d_t1]);  %Note: unitary 2n-inverse-fft 
        mu_Z = u2d_t2 + z_intf; 
    elseif any(contains(["O", "O", "S"], Pcfg))
        mu_Z = Adata * (u1).' + z_intf;

    end

    %W := Ad * U + N1
    var_W = alpha1 * 1 / M * trAdAdh + nu_N1sim*gamma; 

    Nsam = 1; %Number of strings to average 
    w_sam = mu_Z.' + sqrt(var_W / 2) * (randn(Nsam, M) + 1j * randn(Nsam, M)); %

    EQ_log2pYgU = log2(exp(1)) * mean(wdenoise.log_pygw(y, w_sam, nu_p_vec), "all");
    J = DivQP - EQ_log2pYgU; %Cost that is used to adjust adaptive damping 
end
