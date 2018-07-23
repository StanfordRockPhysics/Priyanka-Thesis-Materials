% % load Nishank's computational data
% load('C:\Users\pdutta03\Google Drive\Q-Summer2017\Research\Joint_numerical_data.mat')

% define phase end-points
kf = 2.25; gf = 1e-6; sf = Cw; 
kq = Kqtz; gq = Gqtz; sq = 1e-6; phiq = Phi_sandstone; sqdata = 1.*Cond_sandstone;kqdata = 1.*K_sandstone;
k_mc = Kcal; mu_mc = Gcal; s_mc = 1e-3; phic = Phi_carbonate; sc = Cond_carbonate;

% analyse sandstone samples: conductivity to bulk modulus
% GT conductivity to bulk modulus
[GT1k, GT2k, GT3k, GT4k, GT5k] = runGandT3_PD(sq, sf, kq, kf, gq, gf, phiq, sqdata, 1);
% plot of bulk modulus HS bounds
figure;hold on;xlabel('porosity');ylabel('bulk modulus (GPa)');
[ku,kl,gu,gl,porHS]= hash_PD(kq,gq,kf,gf);
plot(porHS,ku,'-g');plot(porHS,kl,'-g');
% plot computed bulk modulus
plot(phiq,kqdata,'*');
% plot GT predictions
plot(phiq,GT1k,'o');plot(phiq,GT2k,'o');plot(phiq,GT3k,'o');plot(phiq,GT4k,'o');plot(phiq,GT5k,'o');
% plot difference between GT3 and computed K, w.r.t line y = -1
plot(phiq,(GT3k-kqdata'-1),'.');
xL = get(gca,'XLim');
line(xL,[-1 -1]);
% plot self consistent points with aspect ratio 1
[kscq,gscq,por1]=berrysc(kq,gq,kf,gf,1,1);
plot(por1,kscq);
% set transparent legend
legend; set(legend,'color','none');set(legend, 'Box', 'off');
hold off;

% analyse sandstone samples: bulk modulus to conductivity 
% GT bulk modulus to conductivity
[GT1s, GT2s, GT3s, GT4s, GT5s] = runGandT3_PD(sq, sf, kq, kf, gq, gf, phiq, kqdata, 2);
% plot of conductivity HS bounds
figure;hold on;xlabel('porosity');ylabel('conductivity (S/m)');
[su,sl,porHS]= hashc_PD(sq,sf);
plot(porHS,su,'-g');plot(porHS,sl,'-g');
% plot computed conductivity
plot(phiq,sqdata,'*');
% plot GT predictions
plot(phiq,GT1s,'o');plot(phiq,GT2s,'o');plot(phiq,GT3s,'o');plot(phiq,GT4s,'o');plot(phiq,GT5s,'o');
% plot difference between GT3 and computed s, w.r.t line y = -.05
plot(phiq,(GT3s-sqdata'-0.03),'.');
xL = get(gca,'XLim');
line(xL,[-.03 -.03]);
% plot self consistent points with aspect ratio 1
[fscq,por2]=berrysc_res(1/sq,1/sf,1,1,1,1); % fsc = formation factor = Rt/Rw = (in this case) sf/sqdata;
plot(por2,sf./fscq);
% set transparent legend
legend; set(legend,'color','none');set(legend, 'Box', 'off');
hold off;

% plots on bulk modulus - conductivity space
figure;hold on;
plot(su,kl,'-g');plot(sl,ku,'-g');
plot(sqdata,kqdata,'*');
plot(GT1s,kqdata,'o');plot(GT2s,kqdata,'o');plot(GT3s,kqdata,'o');plot(GT4s,kqdata,'o');plot(GT5s,kqdata,'o');
sscq = fliplr(sf./fscq);
plot(sscq, kscq);
xlabel('conductivity');ylabel('bulk modulus (GPa)');
title ('difference between conductivity predicted from bulk modulus and digital conductivity');
plot([sqdata GT3s']',[kqdata kqdata]');
 
 % plots on bulk modulus - conductivity space
figure;hold on;
plot(su,kl,'-g');plot(sl,ku,'-g');
plot(sqdata,kqdata,'*');
plot(sqdata,GT1k,'o');plot(sqdata,GT2k,'o');plot(sqdata,GT3k,'o');plot(GT4s,GT4k,'o');plot(GT5s,GT5k,'o');
sscq = fliplr(sf./fscq);
plot(sscq, kscq);
xlabel('conductivity');ylabel('bulk modulus (GPa)');
title ('difference between bulk modulus predicted from conductivity and digital bulk modulus');
plot([sqdata sqdata]',[kqdata GT3k']');

% select a sandstone point to perform sensitivity analysis (sample 3)
testkq = 19.3731; testsq = 0.0109; testphiq = 0.2083;
% test mineral conductivty end point
sqtemp = [.001+sq sq .002+sq];
figure;
for i = 1:length(sqtemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sqtemp(i), sf, kq, kf, gq, gf, testphiq, testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sqtemp(i), sf, kq, kf, gq, gf, testphiq, testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% test fluid conductivty end point
sftemp = [.9*sf sf 1.1*sf];
figure;title ('fluid conductivity');
for i = 1:length(sftemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sq, sftemp(i), kq, kf, gq, gf, testphiq, testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sq, sftemp(i), kq, kf, gq, gf, testphiq, testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% test mineral bulk end point
kqtemp = [.9*kq kq 1.1*kq];
figure;title ('min bulk');
for i = 1:length(kqtemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sq, sf, kqtemp(i), kf, gq, gf, testphiq, testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sq, sf, kqtemp(i), kf, gq, gf, testphiq, testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% test fluid bulk end point
kftemp = [.9*kf kf 1.1*kf];
figure;
for i = 1:length(kqtemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sq, sf, kq, kftemp(i), gq, gf, testphiq, testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sq, sf, kq, kftemp(i), gq, gf, testphiq, testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% test mineral shear end point
gqtemp = [.9*gq gq 1.1*gq];
figure;
for i = 1:length(kqtemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sq, sf, kq, kf, gqtemp(i), gf, testphiq, testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sq, sf, kq, kf, gqtemp(i), gf, testphiq, testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% test fluid shear end point *************************
gqtemp = [.9*gq gq 1.1*gq];
figure;
for i = 1:length(kqtemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sq, sf, kq, kf, gqtemp(i), gf, testphiq, testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sq, sf, kq, kf, gqtemp(i), gf, testphiq, testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% test porosity end point
testphiqtemp = [.9*testphiq testphiq 1.1*testphiq];
figure;
for i = 1:length(kqtemp)
    [GT1s(i), GT2s(i), GT3s(i), GT4s(i), GT5s(i)] = runGandT3_PD(sq, sf, kq, kf, gq, gf, testphiqtemp(i), testkq, 2);
    [GT1k(i), GT2k(i), GT3k(i), GT4k(i), GT5k(i)] = runGandT3_PD(sq, sf, kq, kf, gq, gf, testphiqtemp(i), testsq, 1);
    subplot(1,3,1);hold on; plot(i, GT1k(i),'o');plot(i, GT2k(i),'o');plot(i, GT3k(i),'o');plot(i, GT4k(i),'o');plot(i, GT5k(i),'o');plot(i, testkq,'*');ylabel('bulk modulus');hold off;
    subplot(1,3,2);hold on; plot(i, GT1s(i),'o');plot(i, GT2s(i),'o');plot(i, GT3s(i),'o');plot(i, GT4s(i),'o');plot(i, GT5s(i),'o');plot(i, testsq,'*');ylabel('conductivity');hold off;
    subplot(1,3,3);hold on; plot(GT1s(i),GT1k(i),'o');plot(GT2s(i),GT2k(i),'o');plot(GT3s(i),GT3k(i),'o');plot(GT4s(i),GT4k(i),'o');plot(GT5s(i),GT5k(i),'o');plot(testsq,testkq,'*');xlabel('conductivity');ylabel('bulk modulus');hold off;
end


% % how do points inside the HS bounds plot on the bulk modulus-conductivity
% % space
% % ku;kl,su, sl;
% kfoam = ku; sfoam = sl;
% ksusp = kl; ssusp = su;
% frac = (0:.05:1)';
% kcomp = frac*kfoam + (1-frac)*ksusp;
% scomp = frac*sfoam + (1-frac)*ssusp;
% figure;hold on;
% xlabel('conductivity');ylabel('bulk modulus (GPa)');
% col = frac*(0.*kfoam+1);
% col2 = (0.*frac+1)*porHS;
% scatter(scomp(:),kcomp(:),15,col2(:));
% t1 = scomp(:);t2 = kcomp(:);t3 = col2(:);
% figure;plot(t1(t3==.5),t2(t3==.5),'o');hold on;
% plot(t1(t3==.1),t2(t3==.1),'o');
% plot(t1(t3==.9),t2(t3==.9),'o');
% plot(t1(t3==.7),t2(t3==.7),'o');
% plot(t1(t3==.3),t2(t3==.3),'o');
% plot(sqdata,kqdata,'*');


% scatter plot of K-sigma data colored by porosity
phit = [.05 .1 .17 .2 .24 .3 .35]; % test porosity values based on data
figure; hold on;xlabel('conductivity (S/m)');ylabel('bulk modulus (GPa)');
title('embedded bound lenses and data points colored by porosity');
scatter(sqdata,kqdata,30,phiq,'filled');col = colormap('jet');colorbar;
for j = 1:length(phit)
    
    [HS1_s,HS2_s]= HSaverageCondB(1-phit(j),sq,sf);
    dx = (HS2_s-HS1_s)/100;
    eps=0;
    s = HS2_s:-dx:HS1_s;
    
    for i = 1:length(s)
        [GT1(i),GT2(i),GT3(i),GT4(i),GT5(i)]=GT((1-phit(j)), phit(j), sq, sf, kq, kf, gq, gf, s(i));
    end
    
    x = (phit(j)-min(phiq))/(max(phiq)-min(phiq));
    plot(s,GT2,'color',col(round(x*63+1),:));
    plot(s,GT3,'color',col(round(x*63+1),:));
    X=[s,fliplr(s)];                %#create continuous x value array for plotting
    Y=[GT2,fliplr(GT3)];              %#create y values for out and then back
    fill(X,Y,col(round(x*63+1),:));
    alpha(.5);
end
scatter(sqdata,kqdata,50,phiq,'filled');colormap('jet');
    
function [FF1,FF2,FF3,FF4,FF5]=GT(f1, f2, s1, s2, k1, k2, g1, g2,data)
a1 = 6.*(g1-g2).*((f1.*s2+f2.*s1+2*s2).^2).*((k1-k2).^2)./( ((s1-s2).^3).*((3*f1.*k2+3*f2.*k1+4*g2).^2) );
a2 = a1.*3*s1.*(3*k1+4*g2)./((s1+2*s2).*(3*k1+4*g1));
a3 = a1.*(2*s1+s2).*(3*k2+4*g2)./(3*s2.*(3*k2+4*g1));
a4 = a1.*2*s1.*g2./((s1+s2).*g1);
a5 = a1.*(s1+s2).*g2./(2.*s2.*g1);
FF1 = makeFF(a1,s1,s2,k1,k2,g1,g2,f1,data);
FF2 = makeFF(a2,s1,s2,k1,k2,g1,g2,f1,data);
FF3 = makeFF(a3,s1,s2,k1,k2,g1,g2,f1,data);
FF4 = makeFF(a4,s1,s2,k1,k2,g1,g2,f1,data);
FF5 = makeFF(a5,s1,s2,k1,k2,g1,g2,f1,data);
end

function FF = makeFF(a,s1,s2,k1,k2,g1,g2,f1,s)

[kHS1,kHS2] = HSaverageB(f1,k1,k2,g1,g2);
[sHS1,sHS2] = HSaverageCondB(f1,s1,s2);

choice  = 1;
if choice == 1
    FF = (a.*kHS1.*(sHS2-s).*(sHS1-sHS2) - kHS2.*(sHS1-s).*(kHS1-kHS2) )./ (a.*(sHS2-s).*(sHS1-sHS2) - (sHS1-s).*(kHS1-kHS2) ); % bulk modulus from conductivity
elseif choice == 2
    FF = -(s*(sHS1*(kHS1 - kHS2) - a*sHS2*(sHS1 - sHS2)) - kHS2*sHS1*(kHS1 - kHS2) + a*kHS1*sHS2*(sHS1 - sHS2))/(kHS2*(kHS1 - kHS2) + s*(kHS2 - kHS1 + a*(sHS1 - sHS2)) - a*kHS1*(sHS1 - sHS2)); % conductivity from bulk modulus
end
end    




