function [wlmn] = Wlmn(l,m,n)

% Legendre coefficients

if m>0 || n>0
    wlmn = 0;
else
    fun = @(ksi,psi,phi) W(ksi,psi,phi)*z(l,m,n,ksi)*exp(-1i*m*psi)*exp(-1i*n*phi);
    wlmn = (4*pi*pi)^(-1)*integral3(fun,0,2*pi,0,2*pi,-1,1);
end
end

function [zlmn] = z(l,m,n,ksi)

%normalized Legendre functions

if l == 0
    zlmn = sqrt(1/2);
elseif l == 2
    zlmn = sqrt(5/2)*(3*(ksi)^2 - 1)/2;
elseif l == 4
    zlmn = sqrt(9/2)*(35*(ksi)^4 - 30*(ksi)^2 + 3)/8;
else
    error('Case outside range of hexagonal symmetry');
end
end

