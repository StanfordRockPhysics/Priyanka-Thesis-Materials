function [y] = hyp_B(x1, y1, x2, y2, x3, y3, x)
% [x, y] = hyp(x1, y1, x2, y2, x3, y3)
%
% Constructs the segment AB of a hyperbola that passes through the pints
% A=(x1, y1), B=(x2, y2), and C = (x3, y3).

% paramter gam traces out curve by spanning 0<= gam <= 1
gam = ((x - x2)*(x1 - x3))/((x - x3)*(x1 - x2));
y = gam*y1 + (1-gam)*y2  - gam.*(1-gam)*(y1 - y2)^2./((1-gam)*y1 + gam*y2 - y3);



