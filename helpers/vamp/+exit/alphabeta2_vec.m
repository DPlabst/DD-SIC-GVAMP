function varargout = alphabeta2_vec(hnd_lmmse_u2w2, u1, w_star, nu_U2_vec, nu_U2_mism_vec, nu_W2_vec, nu_W2_mism_vec, Np2, iterInf, lmmse_trig, rGVAMP)
    
    % Run LMMSE denoiser for (vector) inputs "nu_U2_vec" and "nu_W2_vec"

    alpha2_vec = zeros(1, length(nu_U2_vec));
    beta2_vec = zeros(1, length(nu_W2_vec));

    if isempty(nu_U2_mism_vec) && isempty(nu_W2_mism_vec)
        nu_U2_mism_vec = nu_U2_vec; 
        nu_W2_mism_vec = nu_W2_vec; 
    end 

    for kk = 1:length(nu_U2_vec)

        r2 = u1(nu_U2_mism_vec(kk)); %Noisy version (can be mismatched to nu_U2_vec which the module assumes) 
        p2 = w_star + sqrt(nu_W2_mism_vec(kk)) * Np2; %Noisy version (can be mismatched to nu_W2_vec which the module assumes) 

        [u2d_hat, alpha2_vec(kk), w2_hat, beta2_vec(kk)] = hnd_lmmse_u2w2(r2, nu_U2_vec(kk), p2, nu_W2_vec(kk), iterInf, lmmse_trig, rGVAMP, []);
    end

    %% Return different numbr of output depending if the function is used with vector inputs or not
    varargout{1} = alpha2_vec; 
    varargout{2} = beta2_vec; 
    if isscalar(nu_U2_vec)
        varargout{3} = u2d_hat; 
        varargout{4} = w2_hat; 
        varargout{5} = r2; 
        varargout{6} = p2; 
    end 


end
