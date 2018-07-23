function [FF1,FF2,FF3,FF4,FF5]=GT_2(f1, f2, s1, s2, k1, k2, g1, g2,data,choice)

% computing Gibiansky and Torquato (GT) (1996,1998) hyperbolae (FF1 to FF5), to compute bulk modulus from conductivity (choice = 1) or vice versa (choice=2)for a 2-phase, 3D, isotropic composite.
% Widest of the 5 segments are the GT bounds.
% Inputs: f1, s1, k1, g1 : volume fraction, electrical conductivity, bulk and shear modulus of phase1
%         f2, s2, k2, g2 : volume fraction, electrical conductivity, bulk and shear modulus of phase2
%         data : given conductvity for predicting bulk modulus (for choice
%         = 1) or given bulk modulus for predicting conductvity (choice =
%         2)


a1 = 6.*(g1-g2).*((f1.*s2+f2.*s1+2*s2).^2).*((k1-k2).^2)./( ((s1-s2).^3).*((3*f1.*k2+3*f2.*k1+4*g2).^2) );
a2 = a1.*3*s1.*(3*k1+4*g2)./((s1+2*s2).*(3*k1+4*g1));
a3 = a1.*(2*s1+s2).*(3*k2+4*g2)./(3*s2.*(3*k2+4*g1));
a4 = a1.*2*s1.*g2./((s1+s2).*g1);
a5 = a1.*(s1+s2).*g2./(2.*s2.*g1);
FF1 = makeFF(a1,s1,s2,k1,k2,g1,g2,f1,data,choice);
FF2 = makeFF(a2,s1,s2,k1,k2,g1,g2,f1,data,choice);
FF3 = makeFF(a3,s1,s2,k1,k2,g1,g2,f1,data,choice);
FF4 = makeFF(a4,s1,s2,k1,k2,g1,g2,f1,data,choice);
FF5 = makeFF(a5,s1,s2,k1,k2,g1,g2,f1,data,choice);
end

function FF = makeFF(a,s1,s2,k1,k2,g1,g2,f1,s,choice)

[kHS1,kHS2] = HSaverageB(f1,k1,k2,g1,g2);
[sHS1,sHS2] = HSaverageCondB(f1,s1,s2);

if choice == 1
    FF = (a.*kHS1.*(sHS2-s).*(sHS1-sHS2) - kHS2.*(sHS1-s).*(kHS1-kHS2) )./ (a.*(sHS2-s).*(sHS1-sHS2) - (sHS1-s).*(kHS1-kHS2) ); % bulk modulus from conductivity
elseif choice == 2
    FF = -(s*(sHS1*(kHS1 - kHS2) - a*sHS2*(sHS1 - sHS2)) - kHS2*sHS1*(kHS1 - kHS2) + a*kHS1*sHS2*(sHS1 - sHS2))/(kHS2*(kHS1 - kHS2) + s*(kHS2 - kHS1 + a*(sHS1 - sHS2)) - a*kHS1*(sHS1 - sHS2)); % conductivity from bulk modulus
end
end    