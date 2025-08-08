function varargout = alpha1_vec(hnd_denoise_u1, u1, nu_U1_vec, nu_U1_mism_vec, rGVAMP)
    
    % Run input denoiser for (vector) inputs "nu_U1_vec"
    if isempty(nu_U1_mism_vec) %Matched case 
        nu_U1_mism_vec = nu_U1_vec; 
    end 

    alpha1_vec = zeros(1, length(nu_U1_vec));

    for kk = 1:length(nu_U1_vec)
        r1 = u1(nu_U1_mism_vec(kk)); %Could be mismatched to nu_U1_vec which the module uses for denoising 
        [u1hat, alpha1_kk] = hnd_denoise_u1(r1.', nu_U1_vec(kk), rGVAMP); %Variance 
        alpha1_vec(kk) = mean(alpha1_kk);
    end
    
    %% Return different numbr of output depending if the function is used with vector inputs or not
    varargout{1} = alpha1_vec; 
    if isscalar(nu_U1_vec)
        varargout{2} = u1hat; 
        varargout{3} = r1; 
    end 
end
