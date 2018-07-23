% script using GT cross-bounds to predict Archie cementation factor 'm'
% from bulk-modulus and conductivity

% written by Priyanka Dutta, 2017

clear;clc;
phi = 0.3;
c1 = 1e-5;c2 = .2;
k1 = 36.6;g1 = 45;k2 = 2.3;g2 =0;

[HS1_s,HS2_s]= HSaverageCondB(1-phi,c1,c2);
dx = (HS2_s-HS1_s)/100;
s = HS2_s:-dx:HS1_s;

n = 0;
    for j = 1:length(s)
            [GT1(j),GT2(j),GT3(j),GT4(j),GT5(j)]= GT_2((1-phi), phi, c1, c2, k1, k2, g1, g2, s(j),1);
            GT6(j) = (GT3(j)+GT1(j))./2; GT7(j) = .8.*GT3(j)+.2.*GT2(j);
            
            k = GT7(j):(GT6(j)-GT7(j))/10:GT6(j);
            for i = 1:length(k)
                n = n+1;
                so(n) = s(j);
                ko(n) = k(i);
                m(n) = (log(s(j)/c2))/log(phi);
            end                
                
    end

figure;hold on;xlabel('conductivity (S/m)');ylabel('bulk modulus (GPa)');title(['sandstone, phi = ' num2str(phi)]);
scatter(so,ko,50,m,'filled');

%% all grid
clear;clc;
c1 = 1e-5;c2 = .2;
k1 = 36.6;g1 = 45;k2 = 2.3;g2 =0;

sall = min(c1,c2):abs(c1-c2)/200:max(c1,c2);
kall = min(k1,k2):abs(k1-k2)/200:max(k1,k2);

phis = 0:.01:.4;
n = 0;

for i=1:length(phis)
    phi = phis(i);
    [HS1_s,HS2_s]= HSaverageCondB(1-phi,c1,c2);
    s = sall(sall<=HS2_s & sall>=HS1_s);

for j = 1:length(s)
    s(j);
            [GT1(j),GT2(j),GT3(j),GT4(j),GT5(j)]= GT_2((1-phi), phi, c1, c2, k1, k2, g1, g2, s(j),1);
            GT6(j) = (GT3(j)+GT1(j))./2 ;
            GT7(j) = .8.*GT3(j)+.2.*GT2(j);
            
            k = kall(kall<=GT6(j) & kall>=GT7(j));
            for l = 1:length(k)
                n = n+1;
                so(n) = s(j);
                ko(n) = k(l);
                m(n) = (log(s(j)/c2))/log(phi);
            end                    
end

phi;
end

p = 0;
for i = 1:length(sall)
    for j = 1:length(kall)
        clear mn;
        mn = m(so==sall(i) & ko==kall(j));
        if ~isempty(mn)
            p = p+1;
            sp(p)=sall(i);
            kp(p)=kall(j);
            minm(p) = min(mn);
            maxm(p) = max(mn);
        end
    end
end

figure;hold on;xlabel('conductivity');ylabel('bulk modulus');title('minimum m');
scatter(sp,kp,40,minm,'filled');

figure;hold on;xlabel('conductivity');ylabel('bulk modulus');title('maximum m');
scatter(sp,kp,40,maxm,'filled');
        
figure;hold on;xlabel('conductivity');ylabel('bulk modulus');title('range m');
scatter(sp,kp,40,maxm-minm,'filled');
        