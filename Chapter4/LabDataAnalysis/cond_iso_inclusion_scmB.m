function [Condeff]=cond_iso_inclusion_scmB(Cond,asp,fracs)
%condscm - Effective electrical conductivity for multi-component ISOTROPIC composite
% [Condeff]=cond_iso_inclusion_scm(Cond,asp,frac)
%
% Effective ELECTRICAL CONDUCTIVITY for ISOTROPIC composite of spheroid
% inclusions using the SELF CONSISTENT METHOD method (SCM) .
% This version allows multiple spheroid constituents.
%	Cond:      Conductivity of the N constituent phases (vector of length N)
%              Each phase is isotropic and has a scalar conductivity.
%	asp:       Vector of aspect ratios for each of the N constituents
%			  < 1 for oblate spheroids; 
%   fracs:     Volume fraction of each phase. Sum(fracs)=1
%	Condeff:   Scalar Effective conductivity 
%

% Written by Nishank Saxena (Shell) (2016)
%  modified by Gary Mavko
%

global Condglobal aspglobal fracsglobal 
Condglobal=Cond;
aspglobal=asp;
fracsglobal=fracs;

% Validate inputs
if length(Cond) ~= length(asp) | length(Cond) ~= length(fracs),
    'input vectors Cond, asp, and frac must all have the same lengths'
    return;
end;
if sum(fracs)~=1,
    'sum of fractions not equal to 1; Normalizing'
end;

[Condeff,Fval,exitflag] = fminbnd(@(Condeff) deltacond(Condeff),min(Cond),max(Cond));
if exitflag~=1, 'Cannot find a solution', end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X=deltacond(C)
global Condglobal aspglobal fracsglobal 

% loop over each constituent
for i = 1:length(aspglobal)
    [R(i)] = Rfactor_cond_iso(C, Condglobal(i),aspglobal(i));
end
X = abs(sum(fracsglobal.*(Condglobal - C).*R));
