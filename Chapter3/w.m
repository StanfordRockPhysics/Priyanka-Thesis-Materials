function [wlmn] = w(l,m,n,a)

% Legendre coefficients

if m>0 || n>0
    wlmn = 0;
else
    if l == 0
    fun = @(ksi,psi,phi) (1/(8*pi*pi)).*((a.^2)./(((a*sin(acos(ksi))).^2)+(cos(acos(ksi)).^2)).^(3/2)).*sqrt(1/2).*exp(-1i.*m.*psi).*exp(-1i.*n.*phi);
    elseif l == 2
    fun = @(ksi,psi,phi) (1/(8*pi*pi)).*((a.^2)./(((a*sin(acos(ksi))).^2)+(cos(acos(ksi)).^2)).^(3/2)).*(sqrt(5/2).*((3*(ksi).^2 - 1)./2)).*exp(-1i.*m.*psi).*exp(-1i.*n.*phi);
    elseif l == 4
    fun = @(ksi,psi,phi) (1/(8*pi*pi)).*((a.^2)./(((a*sin(acos(ksi))).^2)+(cos(acos(ksi)).^2)).^(3/2)).*(sqrt(9/2).*((35*(ksi).^4 - 30*(ksi).^2 + 3)./8)).*exp(-1i.*m.*psi).*exp(-1i.*n.*phi);
    else
    error('Case outside range of hexagonal symmetry');
    end
    wlmn = triplequad(fun,-1,1,0,2*pi,0,2*pi);
    wlmn = wlmn/(4*pi*pi);
end
end