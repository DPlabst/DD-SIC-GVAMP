function [S_sc, x0, n_tx, x0_ind, P_tx_sym_measured, mean_U_prior, var_U_prior] = ...
        comm_transmitter(P_tx_i, var_n_pre, N_os, S, mod_alph_name, L_sequ, mod_alph, pX)

    S_sc = sqrt(P_tx_i) .* S;
    x0_ind = randsample(mod_alph.Q, L_sequ, 1, pX);
    x0 = S_sc(x0_ind);

    mean_U_prior = sum(pX .* S_sc);
    var_U_prior = sum(abs(S_sc).^2 .* pX) - abs(mean_U_prior)^2;

    %% ---- Measure transmit powers for debugging -----
    P_tx_sym_measured = 10 * log10(sum(abs(x0).^2) / (length(x0)));

    %% ------- Add cs. Complex AWGN at the TX (pre-intensity noise) -----------
    % Effective bandwidth is [-1,1] x Rsym
    if var_n_pre > 0
        n_tx = sqrt(var_n_pre / 2) * (randn(1, N_os * L_sequ) + 1j * randn(1, N_os * L_sequ)); %Add white CSCG
    else
        n_tx = 0; %Do nothing
    end

end
