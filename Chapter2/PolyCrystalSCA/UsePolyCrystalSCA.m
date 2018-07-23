% Script for poly-crystal SCA model of alpha-quartz, along with
% poly-crystal Voigt-Reuss and Hashin-Shtrikman bounds. the 'c' in this
% script is for quartz, can be modified to reflect 'c' of any other mineral
% (any symmetry), in Voigt notation.

% Written by Priyanka Dutta, 2018

% Voigt notation elastic stiffness tensor of alpha-quartz, trigonal symmetry (RPH)
k(5) = 37.8;g(5) = 44.3;
c11 = 86.6; c12 = 6.7; c13 = 12.6; c33 = 106.1; c44 = 57.8; c66 = (c11-c12)/2; c14 = -17.8;
c = [c11 c12 c13 c14 0 0;c12 c11 c13 -c14 0 0;c13 c13 c33 0 0 0;c14 -c14 0 c44 0 0;0 0 0 0 c44 c14;0 0 0 0 c14 c66];

% Brown's poly-crystal HS bounds
[hsbrn(:,:,1),vrh(:,:,1)]=HSBounds(c); 

% Kube's poly-crystal HS bounds (K,G in rows)
ret = bounds(1,1,c(1,1),c(1,2),c(1,3),c(1,4),c(1,5),c(1,6),c(2,2),c(2,3),c(2,4),c(2,5),c(2,6),c(3,3),c(3,4),c(3,5),c(3,6),c(4,4),c(4,5),c(4,6),c(5,5),c(5,6),c(6,6));
hskub(1,:) = [ret(2,2) ret(1,2)];
ret = bounds(2,1,c(1,1),c(1,2),c(1,3),c(1,4),c(1,5),c(1,6),c(2,2),c(2,3),c(2,4),c(2,5),c(2,6),c(3,3),c(3,4),c(3,5),c(3,6),c(4,4),c(4,5),c(4,6),c(5,5),c(5,6),c(6,6));
hskub(2,:) = [ret(2,2) ret(1,2)];

% Poly-crystal Voigt-Reuss-Hill
[kv,kr,kh,gv,gr,gh] = IsotropicPolycrystalVRH(c);

% aspect ratios for SCA
a = exp(-3:.4:3);

% plot framework
figure;hold on;title('a-quartz:trigonal, * are SCA points');
h1 = gcf;set(gca,'xscale','log');
plot(a,a.*0+kv,'--','DisplayName','Voigt');plot(a,a.*0+kr,'--','DisplayName','Reuss');plot(a,a.*0+kh,'--','Linewidth', 2,'DisplayName','Hill');xlabel('aspect ratio');ylabel('bulk modulus(GPa)');
plot(a,a.*0+hsbrn(1,1),'-o','DisplayName','HSbrn+');plot(a,a.*0+hsbrn(2,1),'-o','DisplayName','HSbrn-');
plot(a,a.*0+hskub(1,1),'-s','DisplayName','HSkub+');plot(a,a.*0+hskub(2,1),'-s','DisplayName','HSkub-');
legend('show');

figure;hold on;title('a-quartz:trigonal, * are SCA points');
h2 = gcf;set(gca,'xscale','log');
plot(a,a.*0+gv,'--','DisplayName','Voigt');plot(a,a.*0+gr,'--','DisplayName','Reuss');plot(a,a.*0+gh,'--','Linewidth', 2,'DisplayName','Hill');xlabel('aspect ratio');ylabel('shear modulus(GPa)');
plot(a,a.*0+hsbrn(1,2),'-o','DisplayName','HSbrn+');plot(a,a.*0+hsbrn(2,2),'-o','DisplayName','HSbrn-');
plot(a,a.*0+hskub(1,2),'-s','DisplayName','HSkub+');plot(a,a.*0+hskub(2,2),'-s','DisplayName','HSkub-');
legend('show');

% SCA model and plots
for n = 1:length(a)
    [ ksc(n), gsc(n) ] = IsotropicPolycrystalSCA( c , a(n));
    figure(h1); plot(a(n),ksc(n),'*','Markersize',8, 'Linewidth',2);
    figure(h2); plot(a(n),gsc(n),'*','Markersize',8, 'Linewidth',2);
end