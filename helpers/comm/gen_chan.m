function A = gen_chan(ch, h_comb, N_span, N_os, N, M)

    % Generate Channel Matrix
    % Only circulant channels supported 
    
    if ch == 0 %Fiber
        h_comb_pad = [h_comb, zeros(1, N_os * N - (N_span * N_os) - 1)];
        Aup = gallery('circul', h_comb_pad).'; %circulant
        A = Aup(:, 1:N_os:end);
    else
        disp('Channel not supported.')

    end
