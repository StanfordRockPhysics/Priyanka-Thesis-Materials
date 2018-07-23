function [kAu, kBu, kAl, kBl, kkAu, kkBu, kkAl, kkBl] = MiltonSubBounds_plot(f1, f2, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B)
% [k1, k2, k3, k4] = MiltonSubBounds(f1, f2, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B, kdata)
% Computes the solid sub bounds of Viogradov and Milton, 2005, J. Mechanics
% and Physics of Solids, 53, 1248-1279
% Initial and final compositions are specified, so that both constituents
% can be changed.
% V&M give four curves, the outer two being the bounds.  All four are
% plotted here.  Horizontal axis of plot is initial composite bulk modulus
% and the vertical axis is the range of compositie bulk moduli that are
% possible after substitution.  
% A particular initial bulk modulus can be specified and the four new
% estimates are superimposed on the plot
%
% Inputs:    1, 2  refer to the material phase;  A, B refer to initial and final composition
%    f1, f2                          volume fractions
%    k1A, g1A, k2A, g2A     moduli for the first composition (A)
%    k1B, g1B, k2B, g2B       moduli for the new composition (B)
%
% Examples:   The following gives Milton's Figure 1a.
%      f1=.5; f2=.5; k1A = 12; g1A = 14; k2A = 100; g2A = 10; k1B = 1; g1B  = 0.6; k2B = 8; g2B = 6; kdata = 29;
%      MiltonSubBounds_plot(.5, .5, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B)
% Examples:   The following gives Milton's Figure 1b.
%      f1=.5; f2=.5; k1A = 2; g1A = 2; k2A = 16; g2A = 4; k1B = 1; g1B = 0.6; k2B = 8; g2B = 6; kdata = 5.2;
%      MiltonSubBounds_plot(.5, .5, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B)
% Examples:   The following gives Milton's Figure 2c.
%      f1=.5; f2=.5; k1A = 33; g1A = 20; k2A = 188; g2A = 34; k1B = 12; g1B = 2; k2B = 7; g2B = 3; kdata =70;
%      MiltonSubBounds_plot(.5, .5, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B);   
% Examples:   The following gives Milton's Figure 2a.
%      f1 = 0.5; f2 = 0.5; k1A = 225; g1A = 2; k2A = 33; g2A = 18; k1B = 1; g1B = .6; k2B = 7; g2B = 10; kdata=65;
%      MiltonSubBounds_plot(.5, .5, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B)
% Example:    Pure fluid substitution a la Gibiansky and Torquato
%      f1 = 0.7; f2 = 0.3; k1A = 36; g1A = 45; k2A = 3; g2A = 0; k1B = 36; g1B = 45; k2B = 1; g2B = 0; kdata=12;
%      MiltonSubBounds_plot(.5, .5, k1A, g1A, k2A, g2A, k1B, g1B, k2B, g2B)

% Extra notes:
% 1. The Y-transform is inverse of the HS bound equation.  Put in an
%    effective K, and it gives the (4/3)G that makes the HS bound equation fall
%    on K.
% 2.  Inverse Y-transform is the HS bound equation.  Put in y=(4/3)*G and it will
%   compute the K.  (4/3)*Gmax yields HS+ and (4/3)*Gmin yields HS-


% Hashin-Shtrikman bounds for first composition (A)
% I note these as 1 and 2 instead of upper and lower, where  1,2 refer to
% which constitutent shear that is used to weight the HS computation.
% The HS bounds are the end-points whre solid substitution is known exactly
z1A = (4/3)*g1A;
kHS1A = (f1/(k1A+z1A)  + f2/(k2A+z1A))^(-1) - z1A;
z2A = (4/3)*g2A;
kHS2A = (f1/(k1A+z2A)  + f2/(k2A+z2A))^(-1) - z2A;
% Hashin-Shtrikman bounds for second composition (B)
z1B = (4/3)*g1B;
kHS1B = (f1/(k1B+z1B)  + f2/(k2B+z1B))^(-1) - z1B;
z2B = (4/3)*g2B;
kHS2B = (f1/(k1B+z2B)  + f2/(k2B+z2B))^(-1) - z2B;


% special points that fall on the first set of bounds
yP3A = g1A + (1/3)*(f1*g2A + f2*g1A);
yP3B = g1B + (1/3)*(f1*g2B + f2*g1B);
yP4A = g2A + (1/3)*(f2*g1A + f1*g2A);
yP4B = g2B + (1/3)*(f2*g1B + f1*g2B)'
yP5A = (4/3)*(f2*g1A + f1*g2A)
yP5B = (4/3)*(f2*g1B + f1*g2B)
yP6A = (4/3)*g1A*g2A/(f1*g1A + f2*g2A)
yP6B = (4/3)*g1B*g2B/(f1*g1B + f2*g2B)
kP3A = Ytransforminverse(f1,f2, k1A, k2A, yP3A);
kP3B = Ytransforminverse(f1,f2, k1B, k2B, yP3B);
kP4A = Ytransforminverse(f1,f2, k1A, k2A, yP4A);
kP4B = Ytransforminverse(f1,f2, k1B, k2B, yP4B);
kP5A = Ytransforminverse(f1,f2, k1A, k2A, yP5A);
kP5B = Ytransforminverse(f1,f2, k1B, k2B, yP5B);
kP6A = Ytransforminverse(f1,f2, k1A, k2A, yP6A);
kP6B = Ytransforminverse(f1,f2, k1B, k2B, yP6B);

% Milton's first set of curves are found in the Y-domain
% transform the HS bounds to the Y-domain
y1A = Ytransform(f1,f2, k1A, k2A, kHS1A)
y2A = Ytransform(f1,f2, k1A, k2A, kHS2A)
y1B = Ytransform(f1,f2, k1B, k2B, kHS1B)
y2B = Ytransform(f1,f2, k1B, k2B, kHS2B)
% create vectors in the Y-domain to span between the HS end points.
% interpolate between point y1 and point y2 using a straigt line
yAu  = linspace(y1A, y2A, 100);
yBu  = linspace(y1B, y2B, 100);
% interpolate between 1/y1 and 1/y2 with a straight line
yAl  = 1./ linspace(1/y1A, 1/y2A, 100);
yBl  = 1./ linspace(1/y1B, 1/y2B, 100);
% transform back to the K domain
kAu  = Ytransforminverse(f1,f2, k1A, k2A, yAu);
kBu   = Ytransforminverse(f1,f2, k1B, k2B, yBu);
kAl  = Ytransforminverse(f1,f2, k1A, k2A, yAl);
kBl   = Ytransforminverse(f1,f2, k1B, k2B, yBl);
% Plot the first pair of curves
% figure; 
plot(kAu, kBu, '-k');
hold on; 
plot(kAl, kBl, '-k');
plot([kP3A, kP4A, kP5A, kP6A], [kP3B, kP4B, kP5B, kP6B], 'ok')
xlabel('KA');  ylabel('KB');

% Milton's second set of curves corresponds to doubly-coated spheres
% Compute and plot  the second pair of curves
yyAu = yfunction( f2, k1A, g1A, g2A, linspace(0, 1, 100));
yyBu  = yfunction( f2, k1B, g1B, g2B, linspace(0, 1, 100));
yyAl = yfunction( f1, k2A, g2A, g1A, linspace(0, 1, 100));
yyBl  = yfunction( f1, k2B, g2B, g1B, linspace(0, 1, 100));
kkAu  = Ytransforminverse(f1,f2, k1A, k2A, yyAu);
kkBu   = Ytransforminverse(f1,f2, k1B, k2B, yyBu);
kkAl  = Ytransforminverse(f1,f2, k1A, k2A, yyAl);
kkBl   = Ytransforminverse(f1,f2, k1B, k2B, yyBl);
% plot second set of bound curves
plot(kkAu, kkBu, ':k');
plot(kkAl, kkBl, ':k');





% % Repeat the calculation but for only a single specified original modulus kdata
% % transform kdata to the Y-domain
% ykdata = Ytransform(f1,f2, k1A, k2A, kdata)
% % find the  final y-values that correspond to ykdata, by interpolating
% % between y-end points
% % First bound method:
% ydataBu = y1B + (y2B-y1B).* (ykdata-y1A)./(y2A-y1A);
% ydataBl  = 1./(1/y1B + (1/y2B-1/y1B).* (1/ykdata-1/y1A)./(1/y2A-1/y1A))
% kdataBu   = Ytransforminverse(f1,f2, k1B, k2B, ydataBu);
% kdataBl   = Ytransforminverse(f1,f2, k1B, k2B, ydataBl);
% plot( [kdata, kdata], [kdataBu, kdataBl], 'or'); fillsym
% 
% % second bound method
% eu = yfunctioninverse(f2, k1A, g1A, g2A, ykdata);
% el  = yfunctioninverse(f1, k2A, g2A, g1A, ykdata);
% yydataBu =  yfunction(f2, k1B, g1B, g2B, eu)
% yydataBl =   yfunction(f1, k2B, g2B, g1B, el)
% kkdataBu   = Ytransforminverse(f1,f2, k1B, k2B, yydataBu)
% kkdataBl    = Ytransforminverse(f1,f2, k1B, k2B, yydataBl)
% 
% plot( [kdata, kdata], [kkdataBu, kkdataBl], 'ob'); fillsym
% 
%-----------------------------------------------------
function  y = Ytransform(f1,f2, k1, k2, keff)
% Milton's Y-transform of Bulk modulus
% Inputs:
%   f1, f2:    volume fractions of phases
%   k1, k2:   bulk moduli of phases
%   keff:       effective bulk modulus of composite

y = -( f2.*k1+f1.*k2) - f1.*f2.*((k1-k2).^2)./(keff - f1.*k1 - f2.*k2);

;%-----------------------------------------------------
function  keff = Ytransforminverse(f1,f2, k1, k2, y)
% Inverse of Milton's Y-transform
% Inputs:
%   f1, f2:    volume fractions of phases
%   k1, k2:   bulk moduli of phases
%   y:           Y-transform value

keff = (f1.*k1 + f2.*k2) - f1.*f2.*((k1-k2).^2)./(y + f2.*k1+f1.*k2);

;%-----------------------------------------------------
function  y = yfunction(f2, k1, g1, g2, e)
%
% Computes equation 2.7

% e = linspace(0, 1, 100);
ynum = f2*g1*(3*k1 + 4* g2)*( 1-e) + (3*k1 + 4*g1)*g2*e;
yden  =  f2*(3*k1 + 4* g2)*( 1-e) + (3*k1 + 4*g1)*e;
y = (4/3)*ynum./yden;

;%-----------------------------------------------------
function e = yfunctioninverse(f2, k1, g1, g2, y)
% computes inverse of equation 2.7

A= (4/3)*f2*g1*(3*k1 + 4* g2);
B = (4/3)*(3*k1 + 4*g1)*g2;
C =  f2*(3*k1 + 4* g2);
D =  (3*k1 + 4*g1);
e = (A-C*y)/(A-B-C*y+D*y);
