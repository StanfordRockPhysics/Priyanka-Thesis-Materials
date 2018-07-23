function [ iso, cp ] = TensorIsotropicProjection( t )
%Isotropic projection of input stiffness/compliance tensor in Kelvin
% notation
% Based on paper by: Browaeys and Chevrot, 2004
% Coded by: Priyanka Dutta, 2018

lambdaH = (1/3)*[1 1 1 0 0 0;1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];
lambdaS = (1/3)*[2 -1 -1 0 0 0;-1 2 -1 0 0 0;-1 -1 2 0 0 0;0 0 0 3 0 0;0 0 0 0 3 0;0 0 0 0 0 3];

th = sum(sum(t.*lambdaH));
ts = (1/5) * sum(sum(t.*lambdaS));

iso = th*lambdaH + ts*lambdaS;

x = [t(1,1) t(2,2) t(3,3) sqrt(2)*t(2,3) sqrt(2)*t(1,3), sqrt(2)*t(1,2), t(4,4), t(5,5), t(6,6), sqrt(2)*t(1,4), sqrt(2)*t(2,5),sqrt(2)*t(3,6), sqrt(2)*t(3,4), sqrt(2)*t(1,5), sqrt(2)*t(2,6),sqrt(2)*t(2,4), sqrt(2)*t(3,5), sqrt(2)*t(1,6), sqrt(2)*t(5,6), sqrt(2)*t(4,6), sqrt(2)*t(4,5)]';

m = [3/15  3/15  3/15  sqrt(2)/15  sqrt(2)/15  sqrt(2)/15  2/15  2/15  2/15;
3/15  3/15  3/15  sqrt(2)/15  sqrt(2)/15  sqrt(2)/15  2/15  2/15  2/15;
3/15  3/15  3/15  sqrt(2)/15  sqrt(2)/15  sqrt(2)/15  2/15  2/15  2/15;
sqrt(2)/15  sqrt(2)/15  sqrt(2)/15  4/15  4/15  4/15  -sqrt(2)/15  -sqrt(2)/15  -sqrt(2)/15;
sqrt(2)/15  sqrt(2)/15  sqrt(2)/15  4/15  4/15  4/15  -sqrt(2)/15  -sqrt(2)/15  -sqrt(2)/15;
sqrt(2)/15  sqrt(2)/15  sqrt(2)/15  4/15  4/15  4/15  -sqrt(2)/15  -sqrt(2)/15  -sqrt(2)/15;
2/15  2/15  2/15  -sqrt(2)/15  -sqrt(2)/15  -sqrt(2)/15  1/5  1/5  1/5;
2/15  2/15  2/15  -sqrt(2)/15  -sqrt(2)/15  -sqrt(2)/15  1/5  1/5  1/5;
2/15  2/15  2/15  -sqrt(2)/15  -sqrt(2)/15  -sqrt(2)/15  1/5  1/5  1/5];

p = cat(2,m,zeros(9,12));   
p = cat(1,p,zeros(12,21));

mul = [1 1 1 sqrt(2) sqrt(2) sqrt(2) 1 1 1 sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2)]';

iso2 = p*x;
iso2 = iso2.*(1./mul);
cp = [iso2(1)  iso2(6)  iso2(5)  iso2(10) iso2(14) iso2(18);
      iso2(6)  iso2(2)  iso2(4)  iso2(16) iso2(11) iso2(15);
      iso2(5)  iso2(4)  iso2(3)  iso2(13) iso2(17) iso2(12);
      iso2(10) iso2(16) iso2(13) iso2(7)  iso2(21) iso2(20);
      iso2(14) iso2(11) iso2(17) iso2(21) iso2(8)  iso2(19);
      iso2(18) iso2(15) iso2(12) iso2(20) iso2(19) iso2(9)];
end

