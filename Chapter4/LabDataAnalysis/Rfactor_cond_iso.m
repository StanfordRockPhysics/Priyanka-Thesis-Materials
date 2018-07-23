function [R] = Rfactor_cond_iso(CondB, CondI,asp)
% Computes the electric field concentration tensor, R, for a spheroid with
% isotropic background, isotropic fill, and isotropic distribution of
% spheroids
%
% Inputs: 
%         CondB = Background  anisotropic conductivity in the form [CondB11 CondB22 CondB33]
%         CondI = Inclusion   anisotropic conductivity in the form [CondI11 CondI22 CondI33]
%         asp   = Aspect ratio (c/a); where a1 = a2 = a; a3 = c along axes
%                 x1, x2 and x3 (Note if asp = 0, you will model inclusions as layers;
%                 <1 for oblate spheroids
% Output: R     = Electric field concentration tensor (scalar)

% Written by Gary Mavko

gamma = 1/asp;

if gamma > 1,     % oblate spheroid
    A33 = (gamma^2/(gamma^2-1))*(1 - (atan(sqrt(gamma^2-1)))/sqrt(gamma^2-1));
elseif gamma < 1, % prolate spheroid
    A33 = (gamma^2/(gamma^2-1))*(1 - 0.5*log((1+sqrt(1-gamma^2))/(1-sqrt(1-gamma^2)))/sqrt(1-gamma^2));
elseif gamma == 1,
    A33 = 1/3;
end;
A11 = (1-A33)/2;

R = 2/(1+((CondI-CondB)/CondB)*A11) + 1/(1+((CondI-CondB)/CondB)*A33);
R = R/3;

