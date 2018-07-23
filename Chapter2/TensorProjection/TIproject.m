function [cp] = TIproject(c)
% TIproject projects a matrix less symmetric than TI to TI.

% c = input 6X6 Voigt notation stiffness matrix
% cp = output 6X6 Voigt notation TI stiffness matrix

% coded by Priyanka Dutta, 2018
% Reference: Decomposition of the elastic tensor and geophysical
% applications; Browaeys and Chevrot, 2004
 

%% converting input Voigt tensor into 21x1 array
y = [c(1,1) c(2,2) c(3,3) c(2,3) c(1,3) c(1,2) c(4,4) c(5,5) c(6,6) c(1,4) c(2,5) c(3,6) c(3,4) c(1,5) c(2,6) c(2,4) c(3,5) c(1,6) c(5,6) c(4,6) c(4,5)]';
mul = [1 1 1 sqrt(2) sqrt(2) sqrt(2) 1 1 1 sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2) sqrt(2)]';
x = y.*mul;


%% projection multiplier
mhex = [3/8           3/8      0  0   0 1/(4*sqrt(2))  0   0     1/4 ;
        3/8           3/8      0  0   0 1/(4*sqrt(2))  0   0     1/4 ;
         0             0       1  0   0       0        0   0      0 ;
         0             0       0 1/2 1/2      0        0   0      0 ;
         0             0       0 1/2 1/2      0        0   0      0 ;
   1/(4*sqrt(2)) 1/(4*sqrt(2)) 0  0   0      3/4       0   0 -1/(2*sqrt(2)) ;
         0             0       0  0   0       0       1/2 1/2     0 ;
         0             0       0  0   0       0       1/2 1/2     0 ;
        1/4           1/4      0  0   0 -1/(2*sqrt(2)) 0   0     1/2];
  
phex = cat(2,mhex,zeros(9,12));   
phex = cat(1,phex,zeros(12,21));


%% output vector (21x1) to 6X6 Voigt stiffness matrix
z = phex*x;
z = z.*(1./mul);
cp = [z(1)  z(6)  z(5)  z(10) z(14) z(18);
      z(6)  z(2)  z(4)  z(16) z(11) z(15);
      z(5)  z(4)  z(3)  z(13) z(17) z(12);
      z(10) z(16) z(13) z(7)  z(21) z(20);
      z(14) z(11) z(17) z(21) z(8)  z(19);
      z(18) z(15) z(12) z(20) z(19) z(9)];

end

