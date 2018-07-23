clear;clc;

% script comparing effective elastic parameters for poly-crystal muscovite
% from Voigt and Reuss schemes using discretized summations and approximate
% forms from spherical harmonics

% defining the stiffness/compliance tensor (muscovite, source to be
% checked)
c11 = 178; c12 = 42.4; c13 = 14.5; c33 = 54.9; c44 = 12.2; c66 = .5*(c11-c12);
c = [c11 c12 c13 0 0 0;c12 c11 c13 0 0 0;c13 c13 c33 0 0 0;0 0 0 c44 0 0;0 0 0 0 c44 0;0 0 0 0 0 c66];
e = (c(1,1)-c(3,3))/(2*c(3,3));
g = (c(6,6)-c(4,4))/(2*c(4,4));
d = ((c(1,3)+c(4,4))^2 - (c(3,3)-c(4,4))^2)/(2*c(3,3)*(c(3,3)-c(4,4)));
s = inv(c);

% refer to paper 'Effect of grain-scale alignment on seismic anisotropy and
% reflectivity of shales', Johansen et. al.

% This work specifically models the compaction ODF given in Johansen's
% paper. It reveals systematic differences between Thomsen parameters
% computed using Voigt and Reuss angular averaging techniques.

% compaction factor varied from 1:10. Above 10, the ODF is almost
% completely oriented and needs very fine discretization of the angles to
% give relaible results, resulting in unnecessarily high computation times.

% estimates by discrete summation
an = 1:10;
sum = 0;
for n = 1:length(an)
    tic
    otr_v = 0;sum = 0; otr_r = 0;
    a = an(n);
    % angle discretization varied dynamically based on compaction factor
    % for more efficient code. Discretization is made finer as compaction
    % factor goes up. Check: ODF triple integration should add up to one
    % (sum). Here sum > 0.9916
    dtheta = pi/(a*25)
    for theta = 0:dtheta:pi
        for xi = 0:dtheta:2*pi
            for phi = 0:dtheta:2*pi
                [R] = johansen_rotation(theta,xi,phi);
                [OP_v] = bond_rotation( c, R, 1 );
                [OP_r] = bond_rotation( s, R, 2 );
                w = a*a/(8*pi*pi*(cos(theta)*cos(theta)+a*a*sin(theta)*sin(theta))^1.5);
                sum = sum + w*sin(theta)*dtheta^3;
                otr_v = otr_v + w.*sin(theta).*(dtheta^3).*OP_v;
                otr_r = otr_r + w.*sin(theta).*(dtheta^3).*OP_r;
            end
        end
    end
    sum
    Voigt = otr_v
    Reuss = inv(otr_r)
    Hill = 0.5.*(Voigt+Reuss);
    
    C33_v(n) = Voigt(3,3);
    C33_r(n) = Reuss(3,3);
    C33_h(n) = Hill(3,3);
    
    C44_v(n) = Voigt(4,4);
    C44_r(n) = Reuss(4,4);
    C44_h(n) = Hill(4,4);

    eps_v(n) = (Voigt(1,1)-Voigt(3,3))/(2*Voigt(3,3));
    eps_r(n) = (Reuss(1,1)-Reuss(3,3))/(2*Reuss(3,3));
    eps_h(n) = (Hill(1,1)-Hill(3,3))/(2*Hill(3,3));

    gam_v(n) = (Voigt(6,6)-Voigt(4,4))/(2*Voigt(4,4));
    gam_r(n) = (Reuss(6,6)-Reuss(4,4))/(2*Reuss(4,4));
    gam_h(n) = (Hill(6,6)-Hill(4,4))/(2*Hill(4,4));

    del_v(n) = ((Voigt(1,3)+Voigt(4,4))^2 - (Voigt(3,3)-Voigt(4,4))^2)/(2*Voigt(3,3)*(Voigt(3,3)-Voigt(4,4)));
    del_r(n) = ((Reuss(1,3)+Reuss(4,4))^2 - (Reuss(3,3)-Reuss(4,4))^2)/(2*Reuss(3,3)*(Reuss(3,3)-Reuss(4,4)));
    del_h(n) = ((Hill(1,3)+Hill(4,4))^2 - (Hill(3,3)-Hill(4,4))^2)/(2*Hill(3,3)*(Hill(3,3)-Hill(4,4)));
    toc
end


for n = 1:length(an)
    a = an(n);
    

eta_r = (eps_r - del_r)./(1+2.*del_r);eta_v = (eps_v - del_v)./(1+2.*del_v);eta_h = (eps_h - del_h)./(1+2.*del_h);


% estimates by approximate summations using spherical harmonics
% Approximate Voigt summation using Legendre polynomials
[c33lv, c44lv, epslv, gamlv, dellv, etalv ] = AppC( c, an );
% Approximate Reuss summation using Legendre polynomials
[c33lr, c44lr, epslr, gamlr, dellr, etalr ] = AppS( c, an );


figure;hold on;plot(an,eta_r,'o');plot(an,etalr);plot(an,eta_v,'*');plot(an,etalv);plot(an,eta_h,'s');xlabel 'compaction factor';ylabel 'eta';
figure;hold on;plot(an,del_r,'o');plot(an,dellr);plot(an,del_v,'*');plot(an,dellv);plot(an,del_h,'s');xlabel 'compaction factor';ylabel 'delta';
figure;hold on;plot(an,gam_r,'o');plot(an,gamlr);plot(an,gam_v,'*');plot(an,gamlv);plot(an,gam_h,'s');xlabel 'compaction factor';ylabel 'gamma';
figure;hold on;plot(an,eps_r,'o');plot(an,epslr);plot(an,eps_v,'*');plot(an,epslv);plot(an,eps_h,'s');xlabel 'compaction factor';ylabel 'epsilon';
figure;hold on;plot(an,C44_r,'o');plot(an,c44lr);plot(an,C44_v,'*');plot(an,c44lv);plot(an,C44_h,'s');xlabel 'compaction factor';ylabel 'C44';
figure;hold on;plot(an,C33_r,'o');plot(an,c33lr);plot(an,C33_v,'*');plot(an,c33lv);plot(an,C33_h,'s');xlabel 'compaction factor';ylabel 'C33';