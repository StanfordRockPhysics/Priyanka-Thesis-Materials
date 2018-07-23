function [su,sl,por]=hashc_PD(s1, s2)
% Hashin-Shtrikman upper and lower bound curves for conductivity
%
%	s1,s2:	Conductivities of the two constituents
%	su,sl:  Upper and lower bounds on conductivity
%   por: volume fraction of material 2
%
% With no output arguments HASH plots the bounds as a function of porosity
% or fraction of phase 2 material. 
%

% Written by P. Dutta

%********* HS upper and lower bounds **************
%
    por = 0:0.01:1; por(1)= 1e-7;
    smax = max(s1,s2);smin = min(s1,s2);
    if smax == s1
        fmax = por;
    else 
        fmax = 1-por;
    end
    fmin = 1-fmax;
	su = smax + fmax./((smin-smax).^(-1)+ (1-fmax)./(3*smax));
    sl = smin + fmin./((smax-smin).^(-1)+ (1-fmin)./(3*smin));
	
if nargout==0
plot(por,su,'-g',por,sl,'-g','linewidth',1);
end;

