function [u1_hat, alpha1_hat_vec, log_b_x1_norm] = denoise_u1(pU, Ualph, r1, nu_U1, rGVAMP)

    %% Stable implementation in log domain: COMPLEX /REAL
    arg_b_x1 = -abs(r1 - Ualph.').^2 ./ ( (2^(rGVAMP)*nu_U1) ) + log(pU.'); %Projects onto reals if Ualph is real
    log_b_x1_norm = arg_b_x1 - logsumexp(arg_b_x1);
    u1_hat = sum((Ualph.') .* exp(log_b_x1_norm), 1);
    alpha1_hat_vec = (sum((abs(Ualph.').^2) .* exp(log_b_x1_norm), 1) - abs(u1_hat).^2);
    

end
