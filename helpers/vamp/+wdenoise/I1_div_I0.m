function [y,rho_arg] = I1_div_I0(a,p1,nu_Z1,nuN1)
%Compute I1(x) / I0(x) where Ik is the modified Bessel function of the
%first kind, order k 

nu1_prime = nu_Z1 + nuN1; %Concerns only pre-intensity 
rho_arg     =  2*sqrt(a).*abs(p1)./nu1_prime; %This is pre-intensity noise  
bessel_f_sc_I0 = besseli(0,rho_arg,1); 
bessel_f_sc_I1 = besseli(1,rho_arg,1); 

y = exp(log(bessel_f_sc_I1) - log(bessel_f_sc_I0)); 

end

