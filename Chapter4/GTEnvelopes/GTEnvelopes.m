function GTEnvelopes(colmap)

% used to create figures stackedGT3 & stackedGT4

%** quartz and brine properties **%
k1 = 36; g1 = 45; k2 = 2.5; g2 = 0; 
% c1 = 1e-5; c2 = 1/.213;
c1 = .0005; c2 = 4.7;
% k1 = 1; g1 = (3*k1*(1-2*0.3))/(2*(1+0.3)); k2 = 20; g2 = (3*k2*(1-2*0.3))/(2*(1+0.3)); 
% % c1 = 1e-5; c2 = 1/.213;
% c1 = 1; c2 = 20;
f2 = 0:0.05:1; f1 = 1-f2;

%** plots **%

%** rectangles and 5 GT bounds at discreet porosities **%
figure;hold on;
colmap = jet;
phic = 1;
phi = f2;
phi = fliplr([.02,.03,.045,.07,.09,.11,.14,.17,.20,.25,.30,.399]);
%phi = .35
for j=1:length(phi)
    indcolr = round(64*phi(j)/phic);
    colr{j}=num2str(colmap(indcolr,:));
    [k_hsu,k_hsl] = k_hs(1-phi(j),k1,g1,phi(j),k2,g2);
    [HS1_s,HS2_s]= HSaverageCondB(1-phi(j),c1,c2);
    rectangle('Position', [HS1_s,k_hsl,abs(HS2_s-HS1_s),abs(k_hsu-k_hsl)]); %, 'Facecolor',colr{j});
end;

for j=1:length(phi)
%for j=2:length(f2)-1
    % phi = f2;
    %phi = .35
    %indcolr = round(64*phi(j)/phic)
    %colr{j}=num2str(colmap(indcolr,:))
    %[k_hsu,k_hsl] = k_hs(1-phi(j),k1,g1,phi(j),k2,g2);
    if phi == 1
%         phi = phi - 1e-3;
break;
    end
    [HS1_s,HS2_s]= HSaverageCondB(1-phi(j),c1,c2);
    dx = (HS2_s-HS1_s)/100;
    eps=0;
    s = HS2_s:-dx:HS1_s;
    
    phic = 1;
    phie = phi(j)/phic;
    if phie > 1
        phie == 1;
    end
    km = 1./(phic./k2 + (1-phic)./k1);
    kmu = k1 + phie./((km-k1).^-1 + (1-phie).*(k1+4/3.*max(g1,g2)).^-1);
    for i = 1:length(s)
        [GT1(i),GT2(i),GT3(i),GT4(i),GT5(i)]=GT((1-phi(j)), phi(j), c1, c2, k1, k2, g1, g2, s(i));
%         if GT1(i)>kmu
%             GT(i)=kmu;
%             %GT1(i) = [];GT2(i) = [];GT3(i) = [];GT4(i) = [];GT5(i) = [];
%             %break;
%         end
    end
%     GT1(GT1>kmu) = kmu;GT2(GT2>kmu) = kmu;GT3(GT3>kmu) = kmu;GT4(GT4>kmu) = kmu;GT5(GT5>kmu) = kmu;
    %length(s);
    %length(GT1);
   
%     plot(s,GT1,'Color',[.3,.75,.93]);
%     plot(s,GT2,'Color',[.87,.49,.0],'LineWidth',2);
%     plot(s,GT3,'Color',[1,.83,.0],'LineWidth',2);
%     plot(s,GT4,'Color',[.49,.18,.56]);
%     plot(s,GT5,'--','Color',[.47,.67,.19]);      
    plot(s,GT1,'--','Color',[.3,.75,.93]);
    plot(s,GT2,'--','Color',[.87,.49,.0],'LineWidth',2);
    plot(s,GT3,'--','Color',[1,.83,.0],'LineWidth',2);
    plot(s,GT4,'--','Color',[.49,.18,.56]);
    plot(s,GT5,'--','Color',[.47,.67,.19]);
   
end;

f2 = 0:0.005:1; f1 = 1-f2;
%** HS bounds **%
%%** bulk modulus **%
[k_hsu,k_hsl] = k_hs(f1,k1,g1,f2,k2,g2);
%%** conductivity **%
[HS1_s,HS2_s]= HSaverageCondB(1-f2,c1,c2);
%%** HS bounds cross-plots **%
% scatter(HS1_s,k_hsu,50,f2,'DisplayName','HS foam');
% scatter(HS2_s,k_hsu,10,f2,'d','DisplayName','HS suspension-foam');
% scatter(HS1_s,k_hsl,10,f2,'*','DisplayName','HS foam-suspension');
% scatter(HS2_s,k_hsl,50,f2,'+','DisplayName','HS suspension');
% legend('show');
plot(HS1_s,k_hsu,'DisplayName','foam - foam');
plot(HS2_s,k_hsu,'d','DisplayName','suspension(C) - foam(K)');
plot(HS1_s,k_hsl,'*','DisplayName','foam(C) - suspension(K)');
plot(HS2_s,k_hsl,'+','DisplayName','suspension - suspension');
legend('show');

figure;
plot(f2,k_hsu,'DisplayName','HS+');xlabel 'porosity';ylabel 'bulk modulus (GPa)'; hold on; 
plot(f2,k_hsl,'DisplayName','HS-');
figure;
plot(f2,HS1_s,'DisplayName','HS-');xlabel 'porosity';ylabel 'conductivity (S/m)'; hold on; 
plot(f2,HS2_s,'DisplayName','HS+');
legend('show');


% for j=1:length(f2)-1
%     phi = f2;
%     %indcolr = round(64*phi(j)/phic)
%     %colr{j}=num2str(colmap(indcolr,:))
%     %[k_hsu,k_hsl] = k_hs(1-phi(j),k1,g1,phi(j),k2,g2);
%     [HS1_s,HS2_s]= HSaverageCondB(1-phi(j),c1,c2);
%     dx = (HS2_s-HS1_s)/100;
%     s = HS2_s:-dx:HS1_s;
%     
%     phic = 0.4;
%     if phi(j) > phic
%         phic = phi(j);
%     end
%     phie = phi(j)/phic;
%     km = 1./(phic./k2 + (1-phic)./k1);
%     kmu = k1 + phie./((km-k1).^-1 + (1-phie).*(k1+4/3.*max(g1,g2)).^-1);
%     for i = 1:length(s)
%         [GT1(i),GT2(i),GT3(i),GT4(i),GT5(i)]=GT((1-phi(j)), phi(j), c1, c2, k1, k2, g1, g2, s(i));
% %         if GT1(i)>kmu
% %             GT(i)=kmu;
% %             %GT1(i) = [];GT2(i) = [];GT3(i) = [];GT4(i) = [];GT5(i) = [];
% %             %break;
% %         end
%     end
%     GT1(GT1>kmu) = kmu;GT2(GT2>kmu) = kmu;GT3(GT3>kmu) = kmu;GT4(GT4>kmu) = kmu;GT5(GT5>kmu) = kmu;
%     plot(s(1:i),GT1);
%     plot(s(1:i),GT2);
%     plot(s(1:i),GT3);
%     plot(s(1:i),GT4);
%     plot(s(1:i),GT5); 
%     tmp = s(GT1>=kmu);
%     plot(tmp(1),kmu,'o','MarkerFaceColor','k');
%    
% end;

% [k_hsu,k_hsl] = k_hs(f1,k1,g1,f2,k2,g2);
% %%** conductivity **%
% [HS1_s,HS2_s]= HSaverageCondB(1-f2,c1,c2);
% plot(HS1_s,k_hsu,'o','DisplayName','HS foam');
% %plot(k_hsu,HS2_s);
% %plot(k_hsl,HS1_s);
% plot(HS2_s,k_hsl,'*','DisplayName','HS suspension');
% legend('show');
end

function[k_hsu,k_hsl] = k_hs(f1,k1,g1,f2,k2,g2)
% if f1+f2 ~= 0
%     error('sum of f1 & f2 not equal to 1');
% end

%** Hashin-Shtrikman-Walpole bounds for K, RPH p224 **%
k_hsu = k1 + f2./((k2-k1).^-1 + f1.*(k1+4/3.*max(g1,g2)).^-1);
k_hsl = k1 + f2./((k2-k1).^-1 + f1.*(k1+4/3.*min(g1,g2)).^-1);
end

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

choice  =1;
if choice == 1
    FF = (a.*kHS1.*(sHS2-s).*(sHS1-sHS2) - kHS2.*(sHS1-s).*(kHS1-kHS2) )./ (a.*(sHS2-s).*(sHS1-sHS2) - (sHS1-s).*(kHS1-kHS2) ); % bulk modulus from conductivity
elseif choice == 2
    FF = -(s*(sHS1*(kHS1 - kHS2) - a*sHS2*(sHS1 - sHS2)) - kHS2*sHS1*(kHS1 - kHS2) + a*kHS1*sHS2*(sHS1 - sHS2))/(kHS2*(kHS1 - kHS2) + s*(kHS2 - kHS1 + a*(sHS1 - sHS2)) - a*kHS1*(sHS1 - sHS2)); % conductivity from bulk modulus
end
end