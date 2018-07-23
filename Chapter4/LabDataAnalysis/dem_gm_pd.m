function [k,mu]=dem_gm_pd(k1,mu1,k2o,mu2o,asp,phic,por)
%DEM - Effective elastic moduli using Differential Effective Medium
%      formulation.
%
%[K,MU,POR]=DEM(K1,MU1,K2,MU2,ASP,PHIC)
%
%	K1, MU1:	Bulk and shear moduli of background matrix
%	K2, MU2:	Bulk and shear moduli of inclusions
%	ASP:		Aspect ratio of inclusions
%			<1 for oblate spheroids; >1 for prolate spheroids
%	PHIC:		percolation porosity for modified DEM model
%			=1 for usual DEM
%	K, MU:		Effective bulk and shear moduli
%	POR:		Porosity, fraction of phase 2 (should be less than PHIC).
%			For the modified DEM, where phase two is the
%			critical phase, POR is the actual porosity.

% Written by T. Mukerji

% Modified by Gary Mavko to replace material 2 with lower bound
% [k_u,k2,u_u,mu2]=hash_PD2(k1,mu1,k2,mu2,por);
[k_u,k2,u_u,mu2]=bounds(1,[1-phic,phic],[k1,k2o],[mu1,mu2o]);

% Modified by Priyanka to give results for a specific input porosity and
% to use 'hash_PD2' instead of 'bounds'.
tfinal=por./phic;

global DEMINPT;
DEMINPT=ones(1,6);
DEMINPT(1)=k1; 
DEMINPT(2)=mu1; 
DEMINPT(3)=k2; 
DEMINPT(4)=mu2; 
DEMINPT(5)=asp; 
DEMINPT(6)=phic; 

[tout, yout]=ode45m('demyprime_gm',0.00,tfinal,[k1; mu1],max(k1,k2)*1e-15);
clear DEMINPT;

kv=real(yout(:,1)); muv=real(yout(:,2));
k = kv(end);mu = muv(end);

if por>phic
[k_u,k,u_u,mu]=bounds(1,[1-por,por],[k1,k2o],[mu1,mu2o]);
end

phi=phic*tout;
if nargout==0
   figure; plot(phi,kv,'-g',phi,muv,'--r', 'linewidth', 1);
end;
end

