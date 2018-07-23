function [kv,kr,kh,gv,gr,gh] = IsotropicPolycrystalVRH(c)

% IsotropicPolycrystalVRH:  Computes Voigt, Reuss and Hill estimates of
% bulk (kv, kr, kh) and shear (gv, gr, gh) moduli of isotropic polycrystals

% Inputs
% c: 6X6 matrix: Elastic stiffness tensor of single crystal in Voigt notation

% Outputs
% kv, kr, kh: Voigt, Reuss and Hill estimates of bulk moduli of isotropic polycrystals
% gv, gr, gh: Voigt, Reuss and Hill estimates of shear moduli of isotropic polycrystals

% Reference: The Elastic Behaviour of a Crystalline Aggregate, R. Hill,
% 1952 (Proc. Phys. Soc. London A, 65 (1952), pp. 349-354)

% Dependencies: None

% Coded by Priyanka Dutta, 2017

kv = (1/9)*((c(1,1)+c(2,2)+c(3,3))+2*(c(1,2)+c(2,3)+c(3,1)));
gv = (1/15)*((c(1,1)+c(2,2)+c(3,3))-(c(1,2)+c(2,3)+c(3,1))+3*(c(4,4)+c(5,5)+c(6,6)));

s = inv(c);

kr = 1/((s(1,1)+s(2,2)+s(3,3))+2*(s(1,2)+s(2,3)+s(3,1)));
gr = 15/(4*(s(1,1)+s(2,2)+s(3,3))-4*(s(1,2)+s(2,3)+s(3,1))+3*(s(4,4)+s(5,5)+s(6,6)));

kh = (kv+kr)/2;
gh = (gv+gr)/2;

end

