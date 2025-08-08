function [ser] = comm_calc_ser(u1, X_alph_sc, u_ind, data_ind, ind_sic, S_SIC)
    [~, u_ind_threshold] = comm_calc_hd(X_alph_sc, u1); %Hard decision 
    u_data = u_ind(data_ind);
    u_ind_stage_s = u_data(ind_sic:S_SIC:end); % Take data corresponding to current stage
    ser = mean(u_ind_threshold(ind_sic:S_SIC:end) ~= u_ind_stage_s); %Compute SER
end
