% function [V, D] = principal_axes(A)
%
% Compute the principal axes of symmetric second order tensor A.
% 
% A can have two different formats:
% A: symmetric 3 x 3 matrix representing a second order tensor.
% A: a 6 x 1 vector representing a symmetric second order tensor in matrix notation.
%
% D: a 3 x 1 vector containing the eigenvalues of A in ascending order. 
% V: a 3 x 3 matrix of which the columns are the corresponding unit vectors
%    of a righthanded principal axis system.

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

function [V, D] = principal_axes(A)

%% Convert if necessary matrix notation into a symmetric 3x3 matrix.

if size(A,1) == 6 
  INVSQRT2 = sqrt(0.5);

  M = zeros(3,3);
  M(1,1) = A(1);
  M(2,2) = A(2);
  M(3,3) = A(3);

  M(2,3) = INVSQRT2 * A(4);
  M(3,2) = M(2,3);
  M(1,3) = INVSQRT2 * A(5);
  M(3,1) = M(1,3);
  M(1,2) = INVSQRT2 * A(6);
  M(2,1) = M(1,2);
else
  M = A;
end

%% Compute the eigenvalues and eigenvectors.
[V,D] = eig(M);
D = diag(D);

%% Sort the eigenvalues.
[D, I] = sort(D);
V = V(:,I);

%% Ensure that the principal axis system is righthanded.
if det(V) < 0; %dot(cross(V(:,1),V(:,2)), V(:,3)) < 0
  V(:,1) = -V(:,1);
end

end