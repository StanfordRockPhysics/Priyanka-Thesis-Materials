%% monte-carlo simulations
clear;clc;
% gross count
c = 3*1e5;

% C33: random number between 3 and 150
C33 = (150-3).*rand(c,1) + 3; 
% C33 = 24.2.*ones(c,1);

% C44: random number between 3 and 150
C44 = (150-3).*rand(c,1) + 3;
% C44 = 3.7.*ones(c,1);

% EPS: random number between -0.2 and 3
EPS = (3-(-.2)).*rand(c,1) + (-.2);

% GAM: random number between -0.2 and 3
GAM = (3-(-.2)).*rand(c,1) + (-.2);

% DEL: random number between -0.5 and 1
DEL = (1-(-.5)).*rand(c,1) + (-.5);

% screening for tensor elements;
C66 = (2.*GAM + 1).*C44;
C11 = (2.*EPS + 1).*C33;
C13 = sqrt(2.*DEL.*C33.*(C33-C44) + (C33-C44).^2) - C44;
C12 = C11 - 2.*C66;
check1 = C11>abs(C12);
check2 = (C11+C12).*C33 > 2.*C13.*C13;
check3 = imag(C13)==0;
check4 = C33>C44;
check = check1+check2+check3+check4;

% net samples after check
c11 = C11(check == 4);c12 = C12(check == 4);c33 = C33(check == 4);
c44 = C44(check == 4);c66 = C66(check == 4);c13 = C13(check == 4);
eps = EPS(check == 4);gam = GAM(check == 4);del = DEL(check == 4);
c  = length(c33);

% % histogram of inputs
% figure;
% h(1) = subplot(3,2,1);hist(c33);xlabel('C_{33}');ylabel('no. of samples');
% h(2) = subplot(3,2,2);hist(c44);xlabel('C_{44}');%ylabel('no. of samples);
% h(3) = subplot(3,2,3);hist(eps);xlabel('epsilon');ylabel('no. of samples');
% h(4) = subplot(3,2,4);hist(gam);xlabel('gamma');%ylabel('no. of samples');
% h(5) = subplot(3,2,5);hist(del);xlabel('delta');ylabel('no. of samples');
% pos = get(h,'Position');
% new = mean(cellfun(@(v)v(1),pos(1:2)));
% set(h(5),'Position',[new,pos{end}(2:end)]);
% %subplot(2,3,6);hist(eps);xlabel('epsilon');

% V-R averaging with net samples
for n = 1:length(eps)
    C = [c11(n) c12(n) c13(n) 0 0 0;c12(n) c11(n) c13(n) 0 0 0;c13(n) c13(n) c33(n) 0 0 0;0 0 0 c44(n) 0 0;0 0 0 0 c44(n) 0;0 0 0 0 0 c66(n)];
    c1 = C; s1 = inv(c1);
    kv(n) = (1/9)*((c1(1,1)+c1(2,2)+c1(3,3)) +2*(c1(1,2)+c1(2,3)+c1(3,1)));
    kr(n) = 1/((s1(1,1)+s1(2,2)+s1(3,3))+2*(s1(1,2)+s1(2,3)+s1(3,1)));
    gv(n) = (1/15)*((c1(1,1)+c1(2,2)+c1(3,3))-(c1(1,2)+c1(2,3)+c1(3,1))+3*(c1(4,4)+c1(5,5)+c1(6,6)));
    gr(n) = 15/(4*(s1(1,1)+s1(2,2)+s1(3,3))-4*(s1(1,2)+s1(2,3)+s1(3,1))+3*(s1(4,4)+s1(5,5)+s1(6,6)));
end

dfk = (kv - kr)./(kr);
dfp = (kv - kr +(4/3).*(gv-gr))./(kr + (4/3).*gr);
dfg = (gv - gr)./gr;

check11 = dfk<1;
check22 = dfg<1.5;
check33 = dfp<1;
checko = check11+check22+check33;
c11 = c11(checko == 3);c12 = c12(checko == 3);c33 = c33(checko == 3);
c44 = c44(checko == 3);c66 = c66(checko == 3);c13 = c13(checko == 3);
eps = eps(checko == 3);gam = gam(checko == 3);del = del(checko == 3);
dfk = dfk(checko == 3);dfg = dfg(checko == 3);dfp = dfp(checko == 3);
dfk = dfk';dfg = dfg';dfp = dfp';

% histogram of inputs
figure;
h(1) = subplot(3,2,1);hist(c33);xlabel('C_{33}');ylabel('no. of samples');
h(2) = subplot(3,2,2);hist(c44);xlabel('C_{44}');%ylabel('no. of samples);
h(3) = subplot(3,2,3);hist(eps);xlabel('epsilon');ylabel('no. of samples');
h(4) = subplot(3,2,4);hist(gam);xlabel('gamma');%ylabel('no. of samples');
h(5) = subplot(3,2,5);hist(del);xlabel('delta');ylabel('no. of samples');
pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)]);



%% Basic plotting
clr = [.6 .6 1];
%% dfk
dfk_med_c33 = [prctile(dfk(c33>=0 & c33<15),[10 25 50 75 90])' prctile(dfk(c33>=15 & c33<30),[10 25 50 75 90])' prctile(dfk(c33>=30 & c33<45),[10 25 50 75 90])' prctile(dfk(c33>=45 & c33<60),[10 25 50 75 90])' prctile(dfk(c33>=60 & c33<75),[10 25 50 75 90])' prctile(dfk(c33>=75 & c33<90),[10 25 50 75 90])' prctile(dfk(c33>=90 & c33<105),[10 25 50 75 90])' prctile(dfk(c33>=105 & c33<120),[10 25 50 75 90])' prctile(dfk(c33>=120 & c33<135),[10 25 50 75 90])' prctile(dfk(c33>=135 & c33<150),[10 25 50 75 90])'];
%figure;scatterhist(c33(dfk<=3),dfk(dfk<=3),'MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr';xlabel 'C33';
figure; h(1) = subplot(3,2,1);plot(c33(dfk<=3),dfk(dfk<=3),'.','MarkerSize',1,'Color',clr);ylabel({'(Kv - Kr) / Kr  (composite)'});xlabel ({'C33  (domain)'});hold on;
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c33(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c33(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c33(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c33(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c33(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(c332(dfk2<=3),dfk2(dfk2<=3),'o','MarkerSize',1,'Color','red');

dfk_med_c44 = [prctile(dfk(c44>=0 & c44<15),[10 25 50 75 90])' prctile(dfk(c44>=15 & c44<30),[10 25 50 75 90])' prctile(dfk(c44>=30 & c44<45),[10 25 50 75 90])' prctile(dfk(c44>=45 & c44<60),[10 25 50 75 90])' prctile(dfk(c44>=60 & c44<75),[10 25 50 75 90])' prctile(dfk(c44>=75 & c44<90),[10 25 50 75 90])' prctile(dfk(c44>=90 & c44<105),[10 25 50 75 90])' prctile(dfk(c44>=105 & c44<120),[10 25 50 75 90])' prctile(dfk(c44>=120 & c44<135),[10 25 50 75 90])' prctile(dfk(c44>=135 & c44<150),[10 25 50 75 90])'];
%figure;scatterhist(c44(dfk<=3),dfk(dfk<=3),'MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr';xlabel 'C44';
h(2) = subplot(3,2,2);plot(c44(dfk<=3),dfk(dfk<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr  (composite)';xlabel 'C44  (domain)';hold on;
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c44(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c44(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c44(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c44(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfk_med_c44(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(c442(dfk2<=3),dfk2(dfk2<=3),'o','MarkerSize',1,'Color','red');

dfk_med_eps = [prctile(dfk(eps>=-0.2 & eps<0),[10 25 50 75 90])' prctile(dfk(eps>=0 & eps<.3),[10 25 50 75 90])' prctile(dfk(eps>=.3 & eps<.6),[10 25 50 75 90])' prctile(dfk(eps>=.6 & eps<.9),[10 25 50 75 90])' prctile(dfk(eps>=.9 & eps<1.2),[10 25 50 75 90])' prctile(dfk(eps>=1.2 & eps<1.5),[10 25 50 75 90])' prctile(dfk(eps>=1.5 & eps<1.8),[10 25 50 75 90])' prctile(dfk(eps>=1.8 & eps<2.1),[10 25 50 75 90])' prctile(dfk(eps>=2.1 & eps<2.4),[10 25 50 75 90])' prctile(dfk(eps>=2.4 & eps<2.7),[10 25 50 75 90])' prctile(dfk(eps>=2.7 & eps<=3),[10 25 50 75 90])'];
%figure;scatterhist(eps(dfk<=3),dfk(dfk<=3),'MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr';xlabel 'epsilon';
h(3) = subplot(3,2,3);plot(eps(dfk<=3),dfk(dfk<=3),'.','MarkerSize',1,'Color',clr,'DisplayName','simulated data');ylabel '(Kv - Kr) / Kr  (composite)';xlabel 'epsilon  (domain)';hold on;
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_eps(1,:),'LineStyle','--','Color','r','LineWidth',1,'DisplayName','10 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_eps(2,:),'LineStyle','--','Color','k','LineWidth',1,'DisplayName','25 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_eps(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o','DisplayName','median');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_eps(4,:),'LineStyle','--','Color','k','LineWidth',1,'DisplayName','75 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_eps(5,:),'LineStyle','--','Color','r','LineWidth',1,'DisplayName','90 prctile');
legend show;legend('Location','NorthWest');
%plot(eps2(dfk2<=3),dfk2(dfk2<=3),'o','MarkerSize',1,'Color','red');

dfk_med_gam = [prctile(dfk(gam>=-0.2 & gam<0),[10 25 50 75 90])' prctile(dfk(gam>=0 & gam<.1*3),[10 25 50 75 90])' prctile(dfk(gam>=.1*3 & gam<.2*3),[10 25 50 75 90])' prctile(dfk(gam>=.2*3 & gam<.3*3),[10 25 50 75 90])' prctile(dfk(gam>=.3*3 & gam<.4*3),[10 25 50 75 90])' prctile(dfk(gam>=.4*3 & gam<.5*3),[10 25 50 75 90])' prctile(dfk(gam>=.5*3 & gam<.6*3),[10 25 50 75 90])' prctile(dfk(gam>=.6*3 & gam<.7*3),[10 25 50 75 90])' prctile(dfk(gam>=.7*3 & gam<.8*3),[10 25 50 75 90])' prctile(dfk(gam>=.8*3 & gam<2.7),[10 25 50 75 90])' prctile(dfk(gam>=2.7 & gam<=3),[10 25 50 75 90])'];
%figure;scatterhist(gam(dfk<=3),dfk(dfk<=3),'MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr';xlabel 'gamma';
h(4) = subplot(3,2,4);plot(gam(dfk<=3),dfk(dfk<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr  (composite)';xlabel 'gamma  (domain)';hold on;
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_gam(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_gam(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_gam(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_gam(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfk_med_gam(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(gam2(dfk2<=3),dfk2(dfk2<=3),'o','MarkerSize',1,'Color','red');


dfk_med_del = [prctile(dfk(del>=-.5 & del<-.3),[10 25 50 75 90])' prctile(dfk(del>=-.3 & del<-.1),[10 25 50 75 90])' prctile(dfk(del>=-.1 & del<.1),[10 25 50 75 90])' prctile(dfk(del>=.1 & del<.3),[10 25 50 75 90])' prctile(dfk(del>=.3 & del<.5),[10 25 50 75 90])' prctile(dfk(del>=.5 & del<.7),[10 25 50 75 90])' prctile(dfk(del>=.7 & del<.9),[10 25 50 75 90])' prctile(dfk(del>=.9 & del<1.1),[10 25 50 75 90])'];
%figure;scatterhist(del(dfk<=3),dfk(dfk<=3),'MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr';xlabel 'delta';
h(5) = subplot(3,2,5);plot(del(dfk<=3),dfk(dfk<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Kv - Kr) / Kr  (composite)';xlabel 'delta  (domain)';hold on;
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfk_med_del(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfk_med_del(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfk_med_del(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfk_med_del(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfk_med_del(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(del2(dfk2<=3),dfk2(dfk2<=3),'o','MarkerSize',1,'Color','red');

pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)]);



%% dfg
figure;

dfg_med_c33 = [prctile(dfg(c33>=0 & c33<15),[10 25 50 75 90])' prctile(dfg(c33>=15 & c33<30),[10 25 50 75 90])' prctile(dfg(c33>=30 & c33<45),[10 25 50 75 90])' prctile(dfg(c33>=45 & c33<60),[10 25 50 75 90])' prctile(dfg(c33>=60 & c33<75),[10 25 50 75 90])' prctile(dfg(c33>=75 & c33<90),[10 25 50 75 90])' prctile(dfg(c33>=90 & c33<105),[10 25 50 75 90])' prctile(dfg(c33>=105 & c33<120),[10 25 50 75 90])' prctile(dfg(c33>=120 & c33<135),[10 25 50 75 90])' prctile(dfg(c33>=135 & c33<150),[10 25 50 75 90])'];
% figure;scatterhist(c33(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'C33';hold on;
h(1) = subplot(3,2,1);plot(c33(dfg<=3),dfg(dfg<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr  (composite)';xlabel 'C33  (domain)';hold on;
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c33(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c33(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c33(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c33(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c33(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(c332(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfg_med_c44 = [prctile(dfg(c44>=0 & c44<15),[10 25 50 75 90])' prctile(dfg(c44>=15 & c44<30),[10 25 50 75 90])' prctile(dfg(c44>=30 & c44<45),[10 25 50 75 90])' prctile(dfg(c44>=45 & c44<60),[10 25 50 75 90])' prctile(dfg(c44>=60 & c44<75),[10 25 50 75 90])' prctile(dfg(c44>=75 & c44<90),[10 25 50 75 90])' prctile(dfg(c44>=90 & c44<105),[10 25 50 75 90])' prctile(dfg(c44>=105 & c44<120),[10 25 50 75 90])' prctile(dfg(c44>=120 & c44<135),[10 25 50 75 90])' prctile(dfg(c44>=135 & c44<150),[10 25 50 75 90])'];
% figure;scatterhist(c44(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'C44';hold on;
h(2) = subplot(3,2,2);plot(c44(dfg<=3),dfg(dfg<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr  (composite)';xlabel 'C44  (domain)';hold on;
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c44(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c44(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c44(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c44(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfg_med_c44(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(c442(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfg_med_eps = [prctile(dfg(eps>=-.2 & eps<0),[10 25 50 75 90])' prctile(dfg(eps>=0 & eps<.1*3),[10 25 50 75 90])' prctile(dfg(eps>=.1*3 & eps<.2*3),[10 25 50 75 90])' prctile(dfg(eps>=.2*3 & eps<.3*3),[10 25 50 75 90])' prctile(dfg(eps>=.3*3 & eps<.4*3),[10 25 50 75 90])' prctile(dfg(eps>=.4*3 & eps<.5*3),[10 25 50 75 90])' prctile(dfg(eps>=.5*3 & eps<.6*3),[10 25 50 75 90])' prctile(dfg(eps>=.6*3 & eps<.7*3),[10 25 50 75 90])' prctile(dfg(eps>=.7*3 & eps<.8*3),[10 25 50 75 90])' prctile(dfg(eps>=.8*3 & eps<.9*3),[10 25 50 75 90])' prctile(dfg(eps>=.9*3 & eps<=1*3),[10 25 50 75 90])'];
% figure;scatterhist(eps(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'epsilon';hold on;
h(3) = subplot(3,2,3);plot(eps(dfg<=3),dfg(dfg<=3),'.','MarkerSize',1,'Color',clr,'DisplayName','simulated data');ylabel '(Gv - Gr) / Gr  (composite)';xlabel 'epsilon  (domain)';hold on;
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_eps(1,:),'LineStyle','--','Color','r','LineWidth',1,'DisplayName','10 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_eps(2,:),'LineStyle','--','Color','k','LineWidth',1,'DisplayName','25 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_eps(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o','DisplayName','median');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_eps(4,:),'LineStyle','--','Color','k','LineWidth',1,'DisplayName','75 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_eps(5,:),'LineStyle','--','Color','r','LineWidth',1,'DisplayName','90 prctile');
legend show;legend('Location','NorthWest');
%plot(eps2(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfg_med_gam = [prctile(dfg(gam>=-.2 & gam<0),[10 25 50 75 90])' prctile(dfg(gam>=0 & gam<.1*3),[10 25 50 75 90])' prctile(dfg(gam>=.1*3 & gam<.2*3),[10 25 50 75 90])' prctile(dfg(gam>=.2*3 & gam<.3*3),[10 25 50 75 90])' prctile(dfg(gam>=.3*3 & gam<.4*3),[10 25 50 75 90])' prctile(dfg(gam>=.4*3 & gam<.5*3),[10 25 50 75 90])' prctile(dfg(gam>=.5*3 & gam<.6*3),[10 25 50 75 90])' prctile(dfg(gam>=.6*3 & gam<.7*3),[10 25 50 75 90])' prctile(dfg(gam>=.7*3 & gam<.8*3),[10 25 50 75 90])' prctile(dfg(gam>=.8*3 & gam<.9*3),[10 25 50 75 90])' prctile(dfg(gam>=.9*3 & gam<=1*3),[10 25 50 75 90])'];
%figure;scatterhist(gam(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'gamma';
h(4) = subplot(3,2,4);plot(gam(dfg<=3),dfg(dfg<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr  (composite)';xlabel 'gamma  (domain)';hold on;
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_gam(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_gam(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_gam(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_gam(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfg_med_gam(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(gam2(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfg_med_del = [prctile(dfg(del>=-.5 & del<-.3),[10 25 50 75 90])' prctile(dfg(del>=-.3 & del<-.1),[10 25 50 75 90])' prctile(dfg(del>=-.1 & del<.1),[10 25 50 75 90])' prctile(dfg(del>=.1 & del<.3),[10 25 50 75 90])' prctile(dfg(del>=.3 & del<.5),[10 25 50 75 90])' prctile(dfg(del>=.5 & del<.7),[10 25 50 75 90])' prctile(dfg(del>=.7 & del<.9),[10 25 50 75 90])' prctile(dfg(del>=.9 & del<1.1),[10 25 50 75 90])'];
%figure;scatterhist(del(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'delta';
h(5) = subplot(3,2,5);plot(del(dfg<=3),dfg(dfg<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr  (composite)';xlabel 'delta  (domain)';hold on;
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfg_med_del(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfg_med_del(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfg_med_del(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfg_med_del(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfg_med_del(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(del2(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)]);


%% dfp
dfp_med_c33 = [prctile(dfp(c33>=0 & c33<15),[10 25 50 75 90])' prctile(dfp(c33>=15 & c33<30),[10 25 50 75 90])' prctile(dfp(c33>=30 & c33<45),[10 25 50 75 90])' prctile(dfp(c33>=45 & c33<60),[10 25 50 75 90])' prctile(dfp(c33>=60 & c33<75),[10 25 50 75 90])' prctile(dfp(c33>=75 & c33<90),[10 25 50 75 90])' prctile(dfp(c33>=90 & c33<105),[10 25 50 75 90])' prctile(dfp(c33>=105 & c33<120),[10 25 50 75 90])' prctile(dfp(c33>=120 & c33<135),[10 25 50 75 90])' prctile(dfp(c33>=135 & c33<150),[10 25 50 75 90])'];
% figure;scatterhist(c33(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'C33';hold on;
figure;h(1) = subplot(3,2,1);plot(c33(dfp<=3),dfp(dfp<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Pv - Pr) / Pr  (composite)';xlabel 'C33  (domain)';hold on;
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c33(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c33(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c33(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c33(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c33(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(c332(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfp_med_c44 = [prctile(dfp(c44>=0 & c44<15),[10 25 50 75 90])' prctile(dfp(c44>=15 & c44<30),[10 25 50 75 90])' prctile(dfp(c44>=30 & c44<45),[10 25 50 75 90])' prctile(dfp(c44>=45 & c44<60),[10 25 50 75 90])' prctile(dfp(c44>=60 & c44<75),[10 25 50 75 90])' prctile(dfp(c44>=75 & c44<90),[10 25 50 75 90])' prctile(dfp(c44>=90 & c44<105),[10 25 50 75 90])' prctile(dfp(c44>=105 & c44<120),[10 25 50 75 90])' prctile(dfp(c44>=120 & c44<135),[10 25 50 75 90])' prctile(dfp(c44>=135 & c44<150),[10 25 50 75 90])'];
% figure;scatterhist(c44(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'C44';hold on;
h(2) = subplot(3,2,2);plot(c44(dfp<=3),dfp(dfp<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Pv - Pr) / Pr  (composite)';xlabel 'C44  (domain)';hold on;
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c44(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c44(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c44(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c44(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([7.5 22.5 37.5 52.5 67.5 82.5 97.5 112.5 127.5 142.5],dfp_med_c44(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(c442(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfp_med_eps = [prctile(dfp(eps>=-.2 & eps<0),[10 25 50 75 90])' prctile(dfp(eps>=0 & eps<.1*3),[10 25 50 75 90])' prctile(dfp(eps>=.1*3 & eps<.2*3),[10 25 50 75 90])' prctile(dfp(eps>=.2*3 & eps<.3*3),[10 25 50 75 90])' prctile(dfp(eps>=.3*3 & eps<.4*3),[10 25 50 75 90])' prctile(dfp(eps>=.4*3 & eps<.5*3),[10 25 50 75 90])' prctile(dfp(eps>=.5*3 & eps<.6*3),[10 25 50 75 90])' prctile(dfp(eps>=.6*3 & eps<.7*3),[10 25 50 75 90])' prctile(dfp(eps>=.7*3 & eps<.8*3),[10 25 50 75 90])' prctile(dfp(eps>=.8*3 & eps<.9*3),[10 25 50 75 90])' prctile(dfp(eps>=.9*3 & eps<=1*3),[10 25 50 75 90])'];
% figure;scatterhist(eps(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'epsilon';hold on;
h(3) = subplot(3,2,3);plot(eps(dfp<=3),dfp(dfp<=3),'.','MarkerSize',1,'Color',clr,'DisplayName','simulated data');ylabel '(Pv - Pr) / Pr  (composite)';xlabel 'epsilon  (domain)';hold on;
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_eps(1,:),'LineStyle','--','Color','r','LineWidth',1,'DisplayName','10 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_eps(2,:),'LineStyle','--','Color','k','LineWidth',1,'DisplayName','25 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_eps(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o','DisplayName','median');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_eps(4,:),'LineStyle','--','Color','k','LineWidth',1,'DisplayName','75 prctile');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_eps(5,:),'LineStyle','--','Color','r','LineWidth',1,'DisplayName','90 prctile');
legend show;legend('Location','NorthWest');

dfp_med_gam = [prctile(dfp(gam>=-.2 & gam<0),[10 25 50 75 90])' prctile(dfp(gam>=0 & gam<.1*3),[10 25 50 75 90])' prctile(dfp(gam>=.1*3 & gam<.2*3),[10 25 50 75 90])' prctile(dfp(gam>=.2*3 & gam<.3*3),[10 25 50 75 90])' prctile(dfp(gam>=.3*3 & gam<.4*3),[10 25 50 75 90])' prctile(dfp(gam>=.4*3 & gam<.5*3),[10 25 50 75 90])' prctile(dfp(gam>=.5*3 & gam<.6*3),[10 25 50 75 90])' prctile(dfp(gam>=.6*3 & gam<.7*3),[10 25 50 75 90])' prctile(dfp(gam>=.7*3 & gam<.8*3),[10 25 50 75 90])' prctile(dfp(gam>=.8*3 & gam<.9*3),[10 25 50 75 90])' prctile(dfp(gam>=.9*3 & gam<=1*3),[10 25 50 75 90])'];
%figure;scatterhist(gam(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'gamma';
h(4) = subplot(3,2,4);plot(gam(dfp<=3),dfp(dfp<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Pv - Pr) / Pr  (composite)';xlabel 'gamma  (domain)';hold on;
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_gam(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_gam(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_gam(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_gam(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-0.1 .05*3 .15*3 .25*3 .35*3 .45*3 .55*3 .65*3 .75*3 .85*3 .95*3],dfp_med_gam(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(gam2(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

dfp_med_del = [prctile(dfp(del>=-.5 & del<-.3),[10 25 50 75 90])' prctile(dfp(del>=-.3 & del<-.1),[10 25 50 75 90])' prctile(dfp(del>=-.1 & del<.1),[10 25 50 75 90])' prctile(dfp(del>=.1 & del<.3),[10 25 50 75 90])' prctile(dfp(del>=.3 & del<.5),[10 25 50 75 90])' prctile(dfp(del>=.5 & del<.7),[10 25 50 75 90])' prctile(dfp(del>=.7 & del<.9),[10 25 50 75 90])' prctile(dfp(del>=.9 & del<1.1),[10 25 50 75 90])'];
%figure;scatterhist(del(dfg<=3),dfg(dfg<=3),'MarkerSize',1,'Color',clr);ylabel '(Gv - Gr) / Gr';xlabel 'delta';
h(5) = subplot(3,2,5);plot(del(dfp<=3),dfp(dfp<=3),'.','MarkerSize',1,'Color',clr);ylabel '(Pv - Pr) / Pr  (composite)';xlabel 'delta  (domain)';hold on;
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfp_med_del(1,:),'LineStyle','--','Color','r','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfp_med_del(2,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfp_med_del(3,:),'LineStyle','--','Color','k','LineWidth',2,'Marker','o');
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfp_med_del(4,:),'LineStyle','--','Color','k','LineWidth',1);
plot([-.4 -.2 0 .2 .4 .6 .8 .95],dfp_med_del(5,:),'LineStyle','--','Color','r','LineWidth',1);
%plot(del2(dfg2<=3),dfg2(dfg2<=3),'o','MarkerSize',1,'Color','red');

pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(5),'Position',[new,pos{end}(2:end)]);