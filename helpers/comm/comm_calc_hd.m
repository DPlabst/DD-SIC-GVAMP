function [hd,hd_ind] = comm_calc_hd(values, x)
    %Hard decision 
    [~, ind] = min(abs(x-values),[],2);
    hd_ind = ind; 
    hd = values(ind);
end