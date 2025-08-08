function [h_comb, x_L] = gen_fiber_response(x,N_os,h_tx_lim,L_fiber,R_sym)
%This fiber response function assumes that the fiber is purely linear
%This function could be extended to the NLSE

CDFilter = []; %Initialize CD Module
CDFilter.L_fiber = L_fiber; %[m]
CDFilter.R_sym   = R_sym; %Bd or Symb/s
CDFilter.D = 17*1E-12/(1E-9*1E3); %[s/ (m*m)]
CDFilter.c_light = 299792458; %[m/s]
CDFilter.lambda  = 1550E-9; %[m]
CDFilter.beta2  = -CDFilter.D *CDFilter.lambda^2/(2*pi*CDFilter.c_light); %[s^2/m]

if L_fiber ~= 0
    N_CD = length(h_tx_lim);
    deltaf = (N_os*CDFilter.R_sym)/N_CD;
    f = ([0:N_CD-1] - floor(N_CD/2))*deltaf;
    H_CD_omega = exp(+1j * (2*pi*f).^2 * CDFilter.beta2/2 * CDFilter.L_fiber);
    h_comb = ifft(ifftshift(H_CD_omega).*fft(h_tx_lim));
    h_cd   = fftshift(ifft(ifftshift(H_CD_omega))); 

elseif L_fiber == 0
    h_comb = h_tx_lim; %Do nothing
    h_cd   = 1; %Do nothing  
end

%Generate x_L:
x_L = conv(x,h_cd,'same'); 

end

