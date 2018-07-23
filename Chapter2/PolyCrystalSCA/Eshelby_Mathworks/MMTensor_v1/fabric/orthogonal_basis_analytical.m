% function [Q,N,R] = orthogonal_basis_analytical(A, limit_size_to_rank)
%
% Compute the orthogonal tensor basis Q for fabric tensor eigenvalues A 
% using the direct analytical method. 
%
% A: A 3x1 vector containing the eigenvalues of the second order fabric tensor.
%    A will be first be normalized so that sum(A) = 1. From these values the 
%    fabric tensor of the third kind is computed.
%
% limit_size_to_rank: if true (default), then the size of the result is 
%    limited to the rank r. This gives a clear indication of the rank and no
%    memory is wasted on (almost) zero coefficients. If false, the return
%    values have r = 9.
%
% Q: 9 x r coefficient matrix of the orthogonal basis. 
% N: 9 x r vector containing the squared norms of the orthogonal basis
% R: r x r upper triangular matrix connecting the orthogonal to Cowin's
%    basis.
% INVARIANTS: contains the invariants I (should be 0), II and III.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2011, Maarten Moesen                                      % 
% All rights reserved.                                                    %
%                                                                         %
% Author(s): Maarten Moesen                                               %
% E-mail: moesen.maarten@gmail.com                                        %
%                                                                         %
% This file is part of MMTensor.                                          %                                                                        %
%                                                                         %
% Redistribution and use in source and binary forms, with or without      %
% modification, are permitted provided that the following conditions are  %
% met:                                                                    %
%    * Redistributions of source code must retain the above copyright     %
%      notice, this list of conditions and the following disclaimer.      %
%    * Redistributions in binary form must reproduce the above copyright  %
%      notice, this list of conditions and the following disclaimer in    %
%      the documentation and/or other materials provided with the         %
%      distribution.                                                      %
%    * Neither the name of the Katholieke Universiteit Leuven nor the     %
%      names of its contributors may be used to endorse or promote        %
%      products derived from this software without specific prior written %
%      permission.                                                        %
%                                                                         %
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     %
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       %
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A %
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL KATHOLIEKE         %
% UNIVERSITEIT LEUVEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED %
% TO, PROCUREMENT OF SUBSTITUTE GOODS ORSERVICES; LOSS OF USE, DATA, OR   %
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  %
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    %
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      %
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Q,N,R,INVARIANTS] = orthogonal_basis_analytical(A, limit_size_to_rank)

if nargin == 1
  limit_size_to_rank = true;
end
  
%% Initialize general constants
tol = 1e-14;
SQRT2 = sqrt(2);

%% Make A a proper orthogonal eigenvalue tensor

% Ensure that A is a 3x1 vector.
A = reshape(A, 3, 1);

% Normalize A so that trace(A) = 1
A = A ./ sum(A);

% Decompose A into the fabric tensor of the third kind D.
D = 15/2*A - 5/2;

% Precompute useful constants: 
D2 = D.^2;
DC = [D(2)*D(3); D(1)*D(3); D(1)*D(2)];

% Now the following invariants should hold: 
% I = 0
% II < 0 or II == 0 if and only if A(1) == A(2) == A(3)
% 4*II^3+27*III^2 = 0 if and only if A(1) == A(2) or A(1) == A(3) or A(2) == A(3)

nb_equalities = (abs(A(1)-A(2)) < tol) + (abs(A(1)-A(3)) < tol) + (abs(A(2)-A(3)) < tol);

if (nb_equalities == 0)
  basis_rank = 9;
elseif (nb_equalities == 3)
  basis_rank = 2;
else
  basis_rank = 5;
end
 
I = sum(D);
II = sum(DC);
II2 = II^2;
II3 = II*II2;
II4 = II2*II2;
II5 = II^5;
II6 = II3*II3;

III = prod(D);
III2 = III^2;
III3 = III*III2;
III4 = III2*III2;

% Avoid a division by zero if II == 0 (isotropy).
% The limit in case of isotropy is zero.
IIIonII = 0;
if (basis_rank > 2)
  IIIonII = III/II;
end

INVARIANTS = [I, II, III];
%fprintf('Invariants: I= %g, II= %g, III= %g\n', I, II, III);

%% Direct analytical formulation
II56III135 = 56*II3+135*III2;
II8III27 = 8*II3-27*III2;
II4III27 = 4*II3+27*III2;

% CAUTION! Here I introduce an artificial error to avoid 0/0 when II=0 and
% III=0. I rely on these terms only being used in denominators
% where the numerator would be zero too according to L'Hopital's rule.
if (basis_rank == 2)
  II8III27 = 1;
  II56III135 = 1;
end

% Compute the basis Q
Z0 = zeros(3,1);
D0 = ones(3,1); 

if limit_size_to_rank
  Q = zeros(9,basis_rank);
  N = zeros(basis_rank,1);
else
  Q = zeros(9,9);
  N = zeros(9,1);
end

Q(1:3,1) = D0;
Q(4:6,1) = SQRT2*D0;
Q(1:3,2) = 4/3*D0;
Q(4:6,2) = -2/3*SQRT2*D0;
Q(7:9,2) = 2*D0;
N(1) = 9;
N(2) = 20;

if ~limit_size_to_rank || basis_rank > 2
Q(1:3, 3:5) = [2*D, 4/3*D, (D2 + 6/7*IIIonII*D + 28/105*II)];
Q(4:6, 3:5) = SQRT2*[-D, 4/3*D, (D2 + 6/7*IIIonII*D + 91/105*II)];
Q(7:9, 4:5) = [ -2*D, (-9/7*IIIonII*D + 2/5*II)];
N(3:5) = [ -12*II; -56/3*II; 112/35*II2+54/7*III*IIIonII ];
end

if ~limit_size_to_rank || basis_rank > 5
  Q(1:3, 6:8) = [(2*D2 + 3*IIIonII*D + 4/3*II*D0), 2*(54*II2*III*D2 + (27*III2*II-8*II4)*D+ 81*III3+48*III*II3)./II56III135, II*((12*II3+81*III2)*D2 + 16*II4+108*II*III2)./(6*II8III27)];
  Q(4:6, 6:8) = SQRT2*[-(D2 + 3/2*IIIonII*D + 2/3*II*D0), (108*II2*III*D2 + (54*III2*II-16*II4)*D - 81*III3+60*III*II3)./II56III135, II*((12*II3+81*III2)*D2 + 4*II4+27*II*III2)./(6*II8III27)];
  Q(7:9, 7:9) = [ -2*((108*III2*II+16*II4)*D + 81*III3+12*III*II3)./II56III135, -((36*III*II3+243*III3)*D + 8*II5+54*II2*III2)./(6*II8III27), -(2*D2 + 3*IIIonII*D + 4/3*II*D0)];
  N(6:9) = [4*II2+27*III*IIIonII; -2*(32*II6 + 108*III2*II3 - 729*III4) / II56III135; II4III27^2*II/(6*II8III27); (8/3*II2+18*III*IIIonII)]; 
end

if nargout <= 2
  return;
end

% Compute the R-matrix that allows computing Cowin 

if limit_size_to_rank
  R = zeros(basis_rank,basis_rank);
else
  R = zeros(9,9);
end

R(1:2, 1:2) = [1, 2/3; 
              0,   1]; 

% It is only meaningful to compute R up to the theoretical rank of the basis. 
% Otherwise you'd be making linear combinations of zero vectors. Moreover it is
% numerically unstable to compute the coefficients of R for lower ranks.
if basis_rank >= 5
  R(1:5, 3:5) = [0,   0,               0;
                 0,   0,           -II/5; 
                 1, 4/3,               0;           
                 0,   1,   -9/14*IIIonII;
                 0,   0,               1];

  if (basis_rank == 9)
    R(1:9, 6:9) = [     -4/3*II,                         0,                             4*II2/9,                    -8/9*II; 
                              0,                   3/5*III,                              II2/15,                    -4/3*II; 
                   -3/2*IIIonII,                   -2*II/3,                                 III,                 -2*IIIonII;        
                              0,                   -2*II/7,                            3*III/14,               -3/2*IIIonII; 
                              0, -108*II2*III/(II56III135), II*(28*II3+675*III2)/(6*II56III135), 35*II4III27/(3*II56III135);
                              1,                         0,                             -2*II/3,                        4/3;
                              0,                         1,                -18*II2*III/II8III27,        -36*II*III/II8III27;
                              0,                         0,                                   1,                       2/II;
                              0,                         0,                                   0,                          1];
  end

end

end