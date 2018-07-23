%% Load appropriate data 

% all samples: measured elastic/electrical (CSV): THan_USH_ElasticElectricalData.csv
%              mineralogy                  (CSV): THan_USH_MineralogyData.csv
%              measured+mineralogy         (MAT): THan_USH_ElasticElectricalMineralogyData.mat 

% selected low-clay samples for examples in thesis: HanSelected.xlsx, full
% combined details in HanSelected2.xlsx


%% Analyzing data (.mat file)

% % Resistivity varying with confining pressure
% figure;plot(R_2Hz_8,'Color',[255,153,153]);hold on;plot(R_2Hz_15,'Color',[255,102,102]);
% plot(R_2Hz_20,'Color',[255,51,51]);plot(R_2Hz_26,'Color',[204,0,51]);plot(R_2Hz_40,'Color',[153,0,51]);plot(R_2Hz_60,'Color',[51,0,51]);
% xlabel('Sample number');ylabel('Resistivity (ohm-m)');title('Resistivity varying with confining pressure');legend('R-2Hz-8MPa','R-2Hz-15MPa','R-2Hz-20MPa','R-2Hz-26MPa','R-2Hz-40MPa','R-2Hz-60MPa');
% 
% % Conductivity varying with confining pressure
% figure;plot(1./R_2Hz_8,'Color',[255,153,153]);hold on;plot(1./R_2Hz_15,'Color',[255,102,102]);
% plot(1./R_2Hz_20,'Color',[255,51,51]);plot(1./R_2Hz_26,'Color',[204,0,51]);plot(1./R_2Hz_40,'Color',[153,0,51]);plot(1./R_2Hz_60,'Color',[51,0,51]);
% xlabel('Sample number');ylabel('Conductivity (mho/m)');title('Conductivity varying with confining pressure');legend('C-2Hz-8MPa','C-2Hz-15MPa','C-2Hz-20MPa','C-2Hz-26MPa','C-2Hz-40MPa','C-2Hz-60MPa');
% 
% % Vp varying with confining pressure
% figure;plot(Vp_8,'Color',[255,153,153]);hold on;plot(Vp_15,'Color',[255,102,102]);
% plot(Vp_20,'Color',[255,51,51]);plot(Vp_26,'Color',[204,0,51]);plot(Vp_40,'Color',[153,0,51]);plot(Vp_60,'Color',[51,0,51]);
% xlabel('Sample number');ylabel('Resistivity (ohm-m)');title('Resistivity varying with confining pressure');legend('R-2Hz-8MPa','R-2Hz-15MPa','R-2Hz-20MPa','R-2Hz-26MPa','R-2Hz-40MPa','R-2Hz-60MPa');
% 
% % Resisitivity vs Vp color code: Porosity
% figure; scatter(Vp_8,R_2Hz_8,20,Porosity,'filled');hold on;
% scatter(Vp_60,R_2Hz_60,20,Porosity,'filled');
% 
% % Conductivity vs Vp color code: Porosity
% figure; scatter(Vp_8,1./R_2Hz_8,20,Porosity,'filled');hold on;
% scatter(Vp_60,1./R_2Hz_60,20,Porosity,'filled');
% 
% % Conductivity vs Vp color code: Clay
% clay = (Illite+Kaolinite+Chlorite+Smectite)./(100-Porosity); % since data reported such that mineral concentrations+ porosity = 100
% figure; scatter(Vp_8,1./R_2Hz_8,20,clay,'filled');hold on;
% scatter(Vp_60,1./R_2Hz_60,20,clay,'filled');

clay = Clay./(100-Porosity);

%% Additional properties
% clay already computed
feld = Feldspar./(100-Porosity);
qtz = Quartz./(100-Porosity); cal = Calcite./(100-Porosity);
phi = Porosity./100;

kq = 36.6; gq = 45; cq = 1e-5; rhoq = 2.65;%cq = 0.002
kcl = 20.9; gcl = 6.85; ccl = 1/50; rhocl = 2.58;
kca = 77; gca = 32; cca = 1e-5; rhoca = 2.71;
kfe = 74.5; gfe = 33.7; cfe = 1e-5; rhofe = 2.63; % average feldspar properties (labradorite) from RPH
kbr = 2.29; gbr = 0; cbr = 1/0.213; rhobr = 1.025; 

rhomin = (clay*rhocl + qtz*rhoq + feld*rhofe + cal*rhoca);
density = phi.*rhobr + (1-phi).*rhomin;

kv = clay.*kcl + qtz.*kq + feld.*kfe + cal.*kca;gv = clay.*gcl + qtz.*gq + feld.*gfe + cal.*gca;
kr = 1./(clay./kcl + qtz./kq + feld./kfe + cal./kca);gr = 1./(clay./gcl + qtz./gq + feld./gfe + cal./gca);
kmin = 0.5.*(kv+kr);gmin = 0.5.*(gv+gr);

cv = clay.*ccl + qtz.*cq + feld.*cfe + cal.*cca;
cr = 1./(clay./ccl + qtz./cq + feld./cfe + cal./cca);
cmin = 0.5.*(cv+cr);

ksat_60 = 1e-9.*(Vp_60.*Vp_60 - (4/3)*Vs_60.*Vs_60).*density.*1000;
ksat_40 = 1e-9.*(Vp_40.*Vp_40 - (4/3)*Vs_40.*Vs_40).*density.*1000;
%ksat_26 = 1e-9.*(Vp_26.*Vp_26 - (4/3)*Vs_26.*Vs_26).*density.*1000;
%ksat_20 = 1e-9.*(Vp_20.*Vp_20 - (4/3)*Vs_20.*Vs_20).*density.*1000;
%ksat_15 = 1e-9.*(Vp_15.*Vp_15 - (4/3)*Vs_15.*Vs_15).*density.*1000;
ksat_8 = 1e-9.*(Vp_8.*Vp_8 - (4/3)*Vs_8.*Vs_8).*density.*1000;


%% 4 samples selected further (3[E6],7[1VSF],13[4SU],16[W165.7]), for individual GT lens plots

ks = kmin([3,7,13,16]);
gs = gmin([3,7,13,16]);
cs = cmin([3,7,13,16]);
phi = Porosity([3,7,13,16])./100;
qs = Quartz([3,7,13,16]);

ks60 = ksat_60([3,7,13,16]);ks8 = ksat_8([3,7,13,16]);
cs60 = 1./R_2Hz_60([3,7,13,16]);cs8 = 1./R_2Hz_8([3,7,13,16]);

a = exp(-5:.1:5);
for i = 2:length(ks)
i
    for n = 1:length(a)
    an = a(n);
    phic = 1;
    %phic = min(an,1/an);
    k_sca(n) = berryscm([ks(i), kbr],[gs(i), gbr],[1, an],[1-phi(i), phi(i)])
    c_sca(n) = cond_iso_inclusion_scmB([cs(i), cbr],[1, an],[1-phi(i), phi(i)])
    k_dem(n) = dem_gm_pd(ks(i),gs(i),kbr,gbr,an,phic,phi(i))
    k_dem2(n) = dem_gm_pd(ks(i),gs(i),kbr,gbr,an,.4,phi(i))
    % dem1(k1,g1,k2,g2,an,phic,phi);
    c_dem(n) = cond_iso_inclusion_dem0(cs(i),cbr,an,1,phic,phi(i))
    c_dem2(n) = cond_iso_inclusion_dem0(cs(i),cbr,an,1,.4,phi(i))
    end

[HS1_s,HS2_s]= HSaverageCondB(1-phi(i),cs(i),cbr);
dx = (HS2_s-HS1_s)/50;
s = HS2_s:-dx:HS1_s;

    for j = 1:length(s)
            [GT1(j),GT2(j),GT3(j),GT4(j),GT5(j)]= GT_2((1-phi(i)), phi(i), cs(i),cbr, ks(i),kbr, gs(i),gbr, s(j),1);
    end

GT6 = (GT3+GT1)./2; GT7 = .8.*GT3+.2.*GT2; GT6(end) = GT7(end);GT6(1) = GT7(1);
    
figure;hold on;xlabel('conductivity (S/m)');ylabel('bulk modulus (GPa)');title(['sandstone vol. frac. = ' num2str(qs(i)) ', phi = ' num2str(phi(i))]);
plot(s,GT1,'--','Color',[.3,.75,.93]);plot(s,GT2,'--','Color',[.87,.49,.0]);plot(s,GT3,'--','Color',[1,.83,.0]);plot(s,GT4,'--','Color',[.49,.18,.56]);plot(s,GT5,'--','Color',[.47,.67,.19]);
plot(s,GT6);hold on;plot(s,GT7);
plot(cs60(i),ks60(i),'*');plot(cs8(i),ks8(i),'o');
scatter(c_sca,k_sca,30,log(a),'filled');scatter(c_dem,k_dem,30,log(a),'filled');
scatter(c_dem2,k_dem2,30,log(a),'filled');
patch([s,fliplr(s)],[GT7,fliplr(GT6)],'b','facealpha',.05);%alpha(.05);

[k_hsu,k_hsl] = k_hs(1-phi(i),ks(i),gs(i),phi(i),kbr,gbr);
plot(s,s.*0+k_hsu);plot(s,s.*0+k_hsl);
k = k_hsl:(k_hsu-k_hsl)./50:k_hsu;
plot(k.*0+HS2_s, k);plot(k.*0+HS1_s, k);

hold off;
end


%% all selected samples, series plot

% all samples
clear GT1 GT2 GT3 GT4 GT5 GT6 GT7
[phi,I] = sort(Porosity./100);
ks = kmin(I);
cs = cmin(I);
gs = gmin(I);
qs = Quartz(I);

ks60 = ksat_60(I);ks8 = ksat_8(I);
cs60 = 1./R_2Hz_60(I);cs8 = 1./R_2Hz_8(I);

for i = 1:length(ks)
    [GT1(i),GT2(i),GT3(i),GT4(i),GT5(i)]= GT_2((1-phi(i)), phi(i), cs(i), cbr, ks(i), kbr, gs(i), gbr, cs60(i),1); % conductivity to bulk modulus
    [GTo1(i),GTo2(i),GTo3(i),GTo4(i),GTo5(i)]= GT_2((1-phi(i)), phi(i), cs(i), cbr, ks(i), kbr, gs(i), gbr, ks60(i),2); % bulk modulus to conductivity
end


GT6 = .4.*GT3+.6.*GT1; GT7 = .7.*GT3+.3.*GT2; %GT6(end) = GT7(end);GT6(1) = GT7(1);
GTo6 = .6.*GTo3+.4.*GTo1; GTo7 = .6.*GTo3+.4.*GTo2; %GTo6(end) = GTo7(end);GTo6(1) = GTo7(1);

i = 1:1:length(ks);
figure;hold on;xlabel('index');ylabel('bulk modulus (GPa)');
plot(i,GT1,'o');plot(i,GT2,'o');plot(i,GT3,'o');plot(i,GT4,'o');plot(i,GT5,'o');

scatter(i,ks60,30,Clay(I),'filled');scatter(i,ks8,30,Clay(I),'filled');
patch([i,fliplr(i)],[GT7,fliplr(GT6)],'b','facealpha',.05);%alpha(.05);
hold off;

figure;hold on;xlabel('index');ylabel('conductivity (S/m)');
plot(i,GTo1,'o');plot(i,GTo2,'o');plot(i,GTo3,'o');plot(i,GTo4,'o');plot(i,GTo5,'o');

scatter(i,cs60,30,Clay(I),'filled');scatter(i,cs8,30,Clay(I),'filled');
patch([i,fliplr(i)],[GTo7,fliplr(GTo6)],'b','facealpha',.05);%alpha(.05);
hold off;





%% local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FF1,FF2,FF3,FF4,FF5]=GT(f1, f2, s1, s2, k1, k2, g1, g2, sdata)

a1 = 6.*(g1-g2).*((f1.*s2+f2.*s1+2*s2).^2).*((k1-k2).^2)./( ((s1-s2).^3).*((3*f1.*k2+3*f2.*k1+4*g2).^2) );
a2 = a1.*3*s1.*(3*k1+4*g2)./((s1+2*s2).*(3*k1+4*g1));
a3 = a1.*(2*s1+s2).*(3*k2+4*g2)./(3*s2.*(3*k2+4*g1));
a4 = a1.*2*s1.*g2./((s1+s2).*g1);
a5 = a1.*(s1+s2).*g2./(2.*s2.*g1);
FF1 = makeFF(a1,s1,s2,k1,k2,g1,g2,f1,sdata);
FF2 = makeFF(a2,s1,s2,k1,k2,g1,g2,f1,sdata);
FF3 = makeFF(a3,s1,s2,k1,k2,g1,g2,f1,sdata);
FF4 = makeFF(a4,s1,s2,k1,k2,g1,g2,f1,sdata);
FF5 = makeFF(a5,s1,s2,k1,k2,g1,g2,f1,sdata);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FF = makeFF(a,s1,s2,k1,k2,g1,g2,f1,s)

[kHS1,kHS2] = HSaverageB(f1,k1,k2,g1,g2);
[sHS1,sHS2] = HSaverageCondB(f1,s1,s2);

FF = (a.*kHS1.*(sHS2-s).*(sHS1-sHS2) - kHS2.*(sHS1-s).*(kHS1-kHS2) )./ (a.*(sHS2-s).*(sHS1-sHS2) - (sHS1-s).*(kHS1-kHS2) );
end