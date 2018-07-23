function [condeff2]=cond_iso_inclusion_dem0(cond1,cond2o,asp,fracs,phic, por)
clear global;
% [condeff, por]=cond_iso_inclusion_dem(cond1,cond2,asp,fracs,phic)
% Effective ELECTRICAL CONDUCTIVITY for ISOTROPIC composite of spheroid
% inclusions using the Differential Effective Medium method (DEM) .
% This version allows multiple spheroid inclusions, but they must all have the same
% fill material.
%
%   cond1      Conductivity of host material
%	cond2      Conductivity of the inclusion material. All inclusions must be the same
%	asp        Vector of Aspect ratios for the N inclusion types.
%              <1 for oblate spheroids
%   fracs      Vector of volume proportion of each phase. Sum(frac) should be 1.
%   phic       Critical porosity, where HS bound will be computed and used
%              as inclusion material
%   por        porosity (< phic)
%	condeff    Scalar Effective conductivity 
%
% Multiple inclusions are handled in the DEM algorithm by adding a portion
% of each type at each increment of porosity.  The vector fracs determines
% the volume fraction of each inclusion type, so that at porosity phi,
% the porosity of each inclusion type is phi*fracs.  Similarly, if an
% increment of porosity dphi is added at each time step, then dphi*fracs is
% the incremet of each inclusion type added.

% Written by Gary Mavko
% Slightly modified by Priyanka, returning condeff2 instead of condeff

% Check validity of inputs
if length(fracs) ~= length(asp) ,
    'input vectors asp and frac must have the same lengths'
    return;
end;
if length(cond2o) ~= 1 ,
    'input cond2 must be a scalar'
    return;
end;
if phic <0 | phic > 1, 
    'critical porosity must be between 0 and 1'
    return;
end;
if cond1<0 | cond2o<0,
    'conductivities must be positive'
    return;
end;

% Compute the HS  bound at the critical porosity.  This becomes the
% inclusion material
cond2 = 1/(phic/(cond2o+2*cond2o) + (1-phic)/(cond1+2*cond2o)) - 2*cond2o;

global DEMINPT;
cond1;
cond2;
asp;
fracs;
phic;

DEMINPT.cond1 = cond1;   % background material conductivity
DEMINPT.cond2 = cond2;   % inclusion material conductivity
DEMINPT.asp   = asp  ;   % aspect ratios
DEMINPT.fracs = fracs;   % volume portions of each inclusion type
DEMINPT.phic  = phic;    % critical porosity


tfinal=por./phic;

% integrate the ODE of conductivity vs. porosity
[phiout, condeff]=ode45m('demyprimeN_cond_iso0',0.00,tfinal,cond1,max(cond1,cond2)*1e-15);
condeff2 = condeff(end);

if por>phic
condeff2 = 1/(por/(cond2o+2*cond2o) + (1-por)/(cond1+2*cond2o)) - 2*cond2o;
end

phi=phic*phiout;
% plot result if there are no outputs specified.
if nargout==0
   figure; plot(phi,condeff,'-g','linewidth', 1);
end;
clear DEMINPT;

