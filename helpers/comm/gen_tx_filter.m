function [h_tx] = gen_tx_filter(ps_filter_str,alpha_rolloff,N_os,N_span)

%% ------------------- Transmit filter -------------------------
if strcmp(ps_filter_str,'RC')  %FD-Raised Cosine (bandlimited)
    h_tx = rcosdesign(alpha_rolloff,N_span,N_os,'normal');
    h_tx = sqrt(N_os)*h_tx; 

elseif strcmp(ps_filter_str,'RRC') %FD-Root-Raised Cosine (bandlimited)
    h_tx = rcosdesign(alpha_rolloff,N_span,N_os,'sqrt');
    h_tx = sqrt(N_os)*h_tx; 
end

%% Pad the filter to a span of N_span: 
LL = (N_os*N_span+1) - length(h_tx); 
h_tx = [zeros(1,LL/2),h_tx,zeros(1,LL/2)]; 

%% Sanity check for unit energy PS filter: 
%  ------ Sanity Check power constraint -----
if abs(norm(h_tx)^2/N_os - 1) > 20*eps  %P_tx_opt definition: norm(h_model_long)^2/N_os != 1
    error('Unit energy pulse constraint violated.');
end


end

