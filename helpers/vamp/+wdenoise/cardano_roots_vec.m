function r = cardano_roots_vec(poly,tol)

%% From: https://proofwiki.org/wiki/Cardano%27s_Formula
%% Compute cardano roots in a vectorized fashion 
N = size(poly,1); 
a = poly(:,1);
b = poly(:,2);
c = poly(:,3);
d = poly(:,4); 

q = (3*a.*c - b.^2) ./ (9*a.^2);
r = (9*a.*b.*c - 27*a.^2.*d - 2*b.^3) ./ (54*a.^3);

D = r.^2 + q.^3; %Discriminant 

idx_d_real = D > 0;  %one root is real 
idx_d_imag = D <= 0; %all roots are real

root_vec = zeros(N,1); 
S = zeros(N,1); 
T = zeros(N,1);  

%% Only one root is real (simple case) 
Dsqrt = sqrt(D(idx_d_real)); 
S(idx_d_real) = nthroot(r(idx_d_real) + Dsqrt, 3);
T(idx_d_real) = nthroot(r(idx_d_real) - Dsqrt, 3);
root_vec(idx_d_real) = (S(idx_d_real) + T(idx_d_real)) - b(idx_d_real)./(3*a(idx_d_real));
        
%% Trigonometric solutions when D <=0
theta = acos(r(idx_d_imag)./sqrt(-(q(idx_d_imag).^3)));
tmp1 = 2 * sqrt(-q(idx_d_imag)); 

%x1 = tmp1 .* cos(theta/3) - b(idx_d_imag)./(3*a(idx_d_imag));
x2 = tmp1 .* cos((theta + 2*pi)/3) - b(idx_d_imag)./(3*a(idx_d_imag)); %Desired root 
%x3 = tmp1 .* cos((theta + 4*pi)/3) - b(idx_d_imag)./(3*a(idx_d_imag));
%xmat = [x1,x2,x3]; 
%idx_val = xmat < tol; 

root_vec(idx_d_imag) = x2; %Stitch together 
r = root_vec; 

end 
