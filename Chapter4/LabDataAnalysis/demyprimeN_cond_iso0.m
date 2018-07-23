function yprime=demyprimeN_cond_iso0(t,y)
% function yprime=demyprimeN_cond_iso(t,y)
% This is only called by function   cond_iso_inclusion_dem.m
%
% Effective ELECTRICAL CONDUCTIVITY for ISOTROPIC composite of spheroid
% inclusions using the Differential Effective Medium method (DEM) .
% This version allows multiple spheroid inclusions, but they must all have the same
% fill material.

% Written by G. Mavko

global DEMINPT;

cond1  = DEMINPT.cond1;   % background material conductivity
cond2  = DEMINPT.cond2;   % inclusion material conductivity
asp    = DEMINPT.asp;     % aspect ratio
fracs  = DEMINPT.fracs;   % volume portion of each inclusion type
phic   = DEMINPT.phic;    % critical porosity

n = length(asp);

% background conductivity, which varies with each increment of inclusion
cond=y; 
  
% compute electric field concentration factors for each of the inclusion
% types
for j=1:n,
    R(j) = Rfactor_cond_iso(cond, cond2, asp(j));
end;

Crhs = sum(fracs.*(cond2-cond).*R);
yprime(1)=Crhs/(1-t) ;
