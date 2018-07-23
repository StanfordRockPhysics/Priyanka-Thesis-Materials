function WaxmanSmitsBussian_plot(data)
% plots WaxmanSmit and Bussian data vs. bounds
% written by Gary Mavko
% modified by Priyanka

% Extract data from inputs
n1      = size(data.data1,1);             % number of samples
Sfluid1 = data.Sfluid1;              % fluid conductivities
Sgrain1 = data.data1(:,4);           % grain conductivities inverted from Bussian model
m1      = data.data1(:,5);           % Archie m inverted from Bussian model
phi1    = data.data1(:,2);           % W-S porosities
Sdata1  = data.data1(:,7:10);        % W-S measured conductivities
Qe1     = data.data1(:,3);           % W-S Q parameter

n2      = size(data.data2,1);             % number of samples
Sfluid2 = data.Sfluid2;              % fluid conductivities
phi2    = data.data2(:,2);           % W-S porosities
Sdata2  = data.data2(:,6:14);        % W-S measured conductivities
Qe2     = data.data2(:,3);           % W-S Q parameter

titstr={['fluid conductivities: ',num2str(Sfluid1,3)]};

% plot measured Conductivity vs. porosity
figure;  plot(phi1,Sdata1(:,1),'or');     % with first fluid
hold on; plot(phi1,Sdata1(:,2),'og');     % with second fluid
hold on; plot(phi1,Sdata1(:,3),'ob');     % with third fluid
hold on; plot(phi1,Sdata1(:,4),'ok');     % with fourth fluid
xlabel('porosity');
ylabel('Effective Conductivity');
title(['Waxman-Smit data: ',titstr])
fillsym
box on
msa8
fsa16

% plot normalized Conductivity vs. porosity
figure;  plot(phi1,Sdata1(:,1)./Sfluid1(1),'or');     % with first fluid
hold on; plot(phi1,Sdata1(:,2)./Sfluid1(2),'og');     % with second fluid
hold on; plot(phi1,Sdata1(:,3)./Sfluid1(3),'ob');     % with third fluid
hold on; plot(phi1,Sdata1(:,4)./Sfluid1(4),'ok');     % with fourth fluid
xlabel('porosity');
ylabel('Effective Conductivity normalized');
title(['Waxman-Smit data: ',titstr])
fillsym
box on
msa8
fsa16

% plot normalized Conductivity vs. porosity
figure;  plot(phi1,Sdata1(:,1)./Sfluid1(1),'or');     % with first fluid
hold on; plot(phi1,Sdata1(:,2)./Sfluid1(2),'og');     % with second fluid
hold on; plot(phi1,Sdata1(:,3)./Sfluid1(3),'ob');     % with third fluid
hold on; plot(phi1,Sdata1(:,4)./Sfluid1(4),'ok');     % with fourth fluid
xlabel('porosity');
ylabel('Effective Conductivity');
title(['Waxman-Smit data Normalized: ',titstr])
fillsym
box on
msa8
fsa16

% plot normalized Conductivity vs. porosity; colored by Sgrain1
figure;  scatterqq(phi1,Sdata1(:,1)./Sfluid1(1),Sgrain1,([.1 1]),jet);     % with first fluid
hold on; scatterqq(phi1,Sdata1(:,2)./Sfluid1(2),Sgrain1,([.1 1]),jet);     % with second fluid
hold on; scatterqq(phi1,Sdata1(:,3)./Sfluid1(3),Sgrain1,([.1 1]),jet);     % with third fluid
hold on; scatterqq(phi1,Sdata1(:,4)./Sfluid1(4),Sgrain1,([.1 1]),jet);     % with fourth fluid
xlabel('porosity');
ylabel('Effective Conductivity Normalized');
title(['Waxman-Smit data Normalized: ',titstr])
fillsym
cbtitle('Sgrain',14)
box on
msa8
fsa16

% plot normalized vs. porosity; color by Q
figure;  scatterqq(phi1,Sdata1(:,1)./Sfluid1(1),Qe1,([.1 1]),jet);     % with first fluid
hold on; scatterqq(phi1,Sdata1(:,2)./Sfluid1(2),Qe1,([.1 1]),jet);     % with second fluid
hold on; scatterqq(phi1,Sdata1(:,3)./Sfluid1(3),Qe1,([.1 1]),jet);     % with third fluid
hold on; scatterqq(phi1,Sdata1(:,4)./Sfluid1(4),Qe1,([.1 1]),jet);     % with fourth fluid
xlabel('porosity');
ylabel('Effective Conductivity Normalized');
title(['Waxman-Smit data Normalized: ',titstr])
fillsym
cbtitle('Qe1',14)
box on
msa8
fsa16

% % % % figure; plot(Qe1,Sgrain1);
% % % figure; scatterqq(Qe1,Sgrain1,phi1,[.05,.3],jet);
% % % xlabel('Qe');ylabel('Sgrain');
% % % box on
% % % cbtitle('Phi',14)
% % % title(['Waxman-Smit data; Bussian Sgrain'])
% % % msa8
% % % fsa16

% % % figure; plot(Qe1,phi1,'ok')
% % % xlabel('Qe')
% % % ylabel('Porosity')
% % % title('Waxman-Smit data')
% % % box on
% % % fsa16

% plot bounds on data with fluids 1,2,3, predicted from measurements with fluid 4 (highest conductivity)
for j=1:n1
    [S1a(j),S2a(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1(j),Sfluid1(4),Sgrain1(j),Sfluid1(1),Sdata1(j,4));
    [S1b(j),S2b(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1(j),Sfluid1(4),Sgrain1(j),Sfluid1(2),Sdata1(j,4));
    [S1c(j),S2c(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1(j),Sfluid1(4),Sgrain1(j),Sfluid1(3),Sdata1(j,4));
end;
figure; hold on;
hline(1) = plot(phi1,Sdata1(:,4),'xk');               % measurements with fluid 4
hline(2) = plot(phi1,S1a,'or','markersize',4); plot(phi1,S2a,'or')   % bounds for fluid 1
hline(3) = plot(phi1,S1b,'og'); plot(phi1,S2b,'og')   % bounds for fluid 2
hline(4) = plot(phi1,S1c,'om'); plot(phi1,S2c,'om')   % bounds for fluid 3
% PD'modification
colormap(jet);
hline(8) = scatter(phi1,Sdata1(:,4).*Sfluid1(1)./Sfluid1(4),30,Qe1,'filled');
hline(9) = scatter(phi1,Sdata1(:,4).*Sfluid1(2)./Sfluid1(4),30,Qe1,'filled');
hline(10) = scatter(phi1,Sdata1(:,4).*Sfluid1(3)./Sfluid1(4),30,Qe1,'filled');
% back to Gary
set(hline(2:4),'markersize',4);
fillsym
%hline(5) = plot(phi1,Sdata1(:,1),'or');             % measurements with fluid 1
hline(5) = scatter(phi1,Sdata1(:,1),30,Qe1,'filled');             % measurements with fluid 1
hline(6) = plot(phi1,Sdata1(:,2),'og');             % measurements with fluid 2
hline(7) = plot(phi1,Sdata1(:,3),'om');             % measurements with fluid 3 
% set(hline(5:7),'linewidth',2);
xlabel('Porosity');
ylabel('Effective Conductivity');
% set(hline(5:7),'markersize',4);
%set(hline(5:7));
legend(hline,'Measured high conductivity', ...
         'Predicted bounds, fluid 1', ...
         'Predicted bounds, fluid 2', ...
         'Predicted bounds, fluid 3', ...         
         'Measured, with fluid 1', ...
         'Measured, with fluid 2', ...
         'Measured, with fluid 3', ...
         'Archie line, fluid 1', ...    % PD's addition
         'Archie line, fluid 2', ...    % PD's addition
         'Archie line, fluid 3');       % PD's addition
title(['W-S data: Predict 1,2,3 from 4; Bussians Sgrain',titstr])
box on
msa6
fsa16

Sgrain1low = 1e-2;
% plot bounds on data with fluids 1,2,3, predicted from measurements with fluid 4 (highest conductivity)
for j=1:n1
    [S1a(j),S2a(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(1),Sdata1(j,4));
    [S1b(j),S2b(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(2),Sdata1(j,4));
    [S1c(j),S2c(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(3),Sdata1(j,4));
end;
figure; hold on;
hline(1) = plot(phi1,Sdata1(:,4),'xk');               % measurements with fluid 4
hline(2) = plot(phi1,S1a,'or','markersize',4); plot(phi1,S2a,'or')   % bounds for fluid 1
hline(3) = plot(phi1,S1b,'og'); plot(phi1,S2b,'og')   % bounds for fluid 2
hline(4) = plot(phi1,S1c,'om'); plot(phi1,S2c,'om')   % bounds for fluid 3
% PD'modification
colormap(jet);
hline(8) = scatter(phi1,Sdata1(:,4).*Sfluid1(1)./Sfluid1(4),30,Qe1,'filled');
hline(9) = scatter(phi1,Sdata1(:,4).*Sfluid1(2)./Sfluid1(4),30,Qe1,'filled');
hline(10) = scatter(phi1,Sdata1(:,4).*Sfluid1(3)./Sfluid1(4),30,Qe1,'filled');
% back to Gary
set(hline(2:4),'markersize',4);
fillsym
%hline(5) = plot(phi1,Sdata1(:,1),'or');             % measurements with fluid 1
hline(5) = scatter(phi1,Sdata1(:,1),30,Qe1,'filled');             % measurements with fluid 1
hline(6) = plot(phi1,Sdata1(:,2),'og');             % measurements with fluid 2
hline(7) = plot(phi1,Sdata1(:,3),'om');             % measurements with fluid 3
% set(hline(5:7),'linewidth',2);
xlabel('Porosity');
ylabel('Effective Conductivity');
% set(hline(5:7),'markersize',4);
legend(hline,'Measured high conductivity', ...
         'Predicted bounds, fluid 1', ...
         'Predicted bounds, fluid 2', ...
         'Predicted bounds, fluid 3', ...
         'Measured, with fluid 1', ...
         'Measured, with fluid 2', ...
         'Measured, with fluid 3',...
         'Archie line, fluid 1', ...    % PD's addition
         'Archie line, fluid 2',...     % PD's addition
         'Archie line, fluid 3');       % PD's addition
title({['W-S data: Predict 1,2,3 from 4; Sgrain=',num2str(Sgrain1low)];titstr{1}})
box on
msa6
fsa16


Sgrain1low = 1e-4;
% ScatterPlot bounds on data with 2, predicted from measurements with fluid 4 (highest conductivity)
for j=1:n1
    [S1a(j),S2a(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(1),Sdata1(j,4));
    [S1b(j),S2b(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(2),Sdata1(j,4));
    [S1c(j),S2c(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(3),Sdata1(j,4));
end;
figure; hold on;
plot(phi1,S1b,'og'); plot(phi1,S2b,'og')   % bounds for fluid 2
set(hline(2:4),'markersize',4);
fillsym
scatterqq(phi1,Sdata1(:,2),Qe1,[.1 1.5],jet);             % measurements with fluid 2
xlabel('Porosity');
ylabel('Effective Conductivity');
% set(hline(5:7),'markersize',4);
box on
msa8
fsa16

Sgrain1low = 1e-4;
% ScatterPlot bounds on data with 2, predicted from measurements with fluid 4 (highest conductivity)
for j=1:n1
    [S1a(j),S2a(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(1),Sdata1(j,4));
    [S1b(j),S2b(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(2),Sdata1(j,4));
    [S1c(j),S2c(j)] = GT_S2S_hyp_B(1-phi1(j),Sgrain1low,Sfluid1(4),Sgrain1low,Sfluid1(3),Sdata1(j,4));
end;
figure; hold on;
plot(phi1,S1b-S2b,'og'); plot(phi1,S2b-S2b,'og')   % bounds for fluid 2
% set(hline(2:4),'markersize',4);
fillsym
scatterqq(phi1,Sdata1(:,2)-S2b',Qe1,[.1 1.5],jet);             % measurements with fluid 2
xlabel('Porosity');
ylabel('Effective Conductivity');
% set(hline(5:7),'markersize',4);
box on
msa8
fsa16

