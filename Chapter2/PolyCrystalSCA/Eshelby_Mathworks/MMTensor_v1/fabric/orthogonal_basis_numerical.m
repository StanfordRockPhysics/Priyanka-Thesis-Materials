% function [Q,N,R] = orthogonal_basis_numerical(A, limit_size_to_rank)
%
% Compute the orthogonal tensor basis Q for fabric tensor eigenvalues A 
% using a numerical method based on QR factorization of Cowin's system.
% This method should give similar results as orthogonal_basis_analytical.
% It is however more expensive.
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
% INVARIANTS: 1x3 vector containing the invariants I (should be 0), II and III.

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

function [Q,N,R,INVARIANTS] = orthogonal_basis_numerical(A, limit_size_to_rank)

%% limit_size_to_rank is true by default.
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

%% Now the following invariants should hold: 
% I = 0
% II < 0 or II == 0 if and only if A(1) == A(2) == A(3)
% 4*II^3+27*III^2 = 0 if and only if A(1) == A(2) or A(1) == A(3) or A(2) == A(3)

I = sum(D);
II = sum(DC);
III = prod(D);

INVARIANTS = [I, II, III];
%fprintf('Invariants: I= %g, II= %g, III= %g\n', I, II, III);


%% Determine the rank of the basis
nb_equalities = (abs(A(1)-A(2)) < tol) + (abs(A(1)-A(3)) < tol) + (abs(A(2)-A(3)) < tol);

if (nb_equalities == 0)
  basis_rank = 9;
elseif (nb_equalities == 3)
  basis_rank = 2;
else
  basis_rank = 5;
end

%% Compute Cowin's system using the simplifying invariant (I=0)
% With respect to the original system:
% The rows are permuted to a more practical and meaningful column order 
% -> Rows 1 7 are always linearly independent
% -> Rows 1 7 2 8 3 are linearly independent two D's are distinct
% -> All rows are only linearly independent if all three D's are distinct
% Permutation: K = K(:, [1 7 2 8 4 3 5 6 9]);

if limit_size_to_rank
  K = zeros(9,basis_rank);  
else
  K = zeros(9,9);
end

K(1:6, 1) = 1;
K(1:3, 2) = 2;
K(7:9, 2) = K(4:6,1);

if ~limit_size_to_rank || basis_rank > 2
  K(1:3, 3) = 2*D;
  K(4:6, 3) = -D;
  K(1:3, 4) = 2*K(1:3,3);
  K(7:9, 4) = K(4:6,3);
  K(1:3, 5) = D2;
  K(4:6, 5) = DC;
end

if ~limit_size_to_rank || basis_rank > 5
  K(1:3, 6) = 2*D2;
  K(4:6, 6) = [D2(2)+D2(3); D2(3)+D2(1); D2(1)+D2(2)];
  K(1:3, 7) = 2*D.*D2;
  K(4:6, 7) = -III;
  K(1:6, 8) = K(1:6,5).^2;
  K(1:3, 9) = 2*K(1:3,6);
  K(7:9, 9) = K(4:6,6);
end

%% Convert Cowin's system to vector notation.
K = K .* repmat([1; 1; 1; SQRT2; SQRT2; SQRT2; 2; 2; 2], 1, size(K,2));

%% QR-based computation of MM's orthogonalized system 
[Q,R] = qr(K,0);
%N = sign(diag(R));
N = diag(R);
Q = Q.*repmat(N', 9, 1);
R = R.*repmat(1./N, 1, size(R,2));
N = N.^2;

% Here it is important to realize that if the system is rank-deficient, 
% the 7 or 4 final norms should be zero and the numbers in the 7 or 4 last 
% rows of R basically are rounding errors. They should be ignored.

end