% function [E, nu] = computeEngineeringConstantsSqrt2(C)
%
% Compute the engineering constants of stiffness tensor C.
% C is a 6 x 6 stiffness tensor in sqrt2 notation.
%
% Return: E: [E_11 E_22 E_33 G_23 G_31 G_12]
% Return: nu: [nu_23 nu_13 nu_12 nu_32 nu_13 nu_21] = [major minor]
% The largest poisson ratios are the most accurate.

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

function [E, nu] = computeEngineeringConstantsSqrt2(C)

% Tolerance
tol = 1e-6;

% Compute the symmetric part of C
C_sym = 0.5*(C+C');
C_asym = 0.5*(C-C');
C_error = max(max(abs(C_asym)));

if (max(max(abs(C_asym))) > tol)
    fprintf('The stiffness matrix is not symmetric! Maximum deviation: %g\n', C_error);
end

E = zeros(1, 6);
nu = zeros(1, 6);

if (rank(C) < 6)
  return;
end

% Compute the shear moduli
for i = 4:6
    E(i) = 0.5*C(i,i);
end

% Compute compliance matrix:
S = inv(C_sym(1:3,1:3)); % Inverse of a symmetric matrix is symmetric

for i = 1:3
    E(i) = 1/S(i,i);
end


% Define nu: [nu_23 nu_13 nu_12 nu_32 nu_31 nu_21]

% Better accuracy for larger poisson ratio:
nuscale = -diag(E(1:3))*S(1:3, 1:3);

nu(1) = nuscale(2,3);
nu(2) = nuscale(1,3);
nu(3) = nuscale(1,2);
nu(4) = nuscale(3,2);
nu(5) = nuscale(3,1);
nu(6) = nuscale(2,1);
