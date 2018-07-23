function [S1,S2] = GT_S2S_hyp_B(f1,s1A,s2A,s1B,s2B, sdata)
% computes fluid conductivity fluid substitution.  Uses the double-coated
% sphere equations that have been found to be possible bounds on the
% substitution for fluids and solids.  
%
% Inputs
%    f1        solid fraction (fraction of phase 1)
%    s1A,s2A   phase conductivities at the initial composition
%    s1B,s2B   phase conductivities at the final composition.
%              s2B can be a vector of fluid properties
%    sdata     "measured" effective conductivity at initial composition
% Output
%    S1,S2     vector of predicted effective conductivity bounds with each
%              new fluid conductivity
%  Initial conductivity must fall between the HS bounds.

% Gary Mavko, April 2017


titstr = {['porosity= ',num2str(1-f1,2)];...
          ['Cmin1= ',num2str(s1A,2)]; ...
          ['Cmin2= ',num2str(s1B,2)]; ...
          ['Cbrine1= ',num2str(s2A,2)]};

[sHS1A,sHS2A] = HSaverageCondB(f1,s1A,s2A);       
if sdata < sHS1A | sdata > sHS2A
    sHS1A
    sHS2A
    'impossible input data'
    return
end;

      
      
% compute initial and final HS bounds
for j=1:length(s2B)
    [sHS1A,sHS2A] = HSaverageCondB(f1,s1A,s2A);       
    [sHS1B,sHS2B] = HSaverageCondB(f1,s1B,s2B(j));
    % parametric equations that generate the bounding hyperbolas
    S1(j) = hyp_B(sHS1A, sHS1B, sHS2A, sHS2B, s1A, s1B,    sdata);
    S2(j) = hyp_B(sHS1A, sHS1B, sHS2A, sHS2B, s2A, s2B(j), sdata);
end;

if nargout == 0,
figure; hold on;
    plot(s2B, S1, s2B, S2,'linewidth',2);

    % text(min(x1),max(y1),titstr,'fontsize',12);
    xlabel('Fluid conductivity'); ylabel('Effective Conductivity');
    title('Fluid substitution of Conductivity')
    legend('mineral coated','fluid coated')
    fsa16
end;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = hyp_B(x1, y1, x2, y2, x3, y3, x)
% [x, y] = hyp(x1, y1, x2, y2, x3, y3)
%
% Constructs the segment AB of a hyperbola that passes through the pints
% A=(x1, y1), B=(x2, y2), and C = (x3, y3).

% paramter gam traces out curve by spanning 0<= gam <= 1
gam = ((x - x2)*(x1 - x3))/((x - x3)*(x1 - x2));
y = gam*y1 + (1-gam)*y2  - gam.*(1-gam)*(y1 - y2)^2./((1-gam)*y1 + gam*y2 - y3);


