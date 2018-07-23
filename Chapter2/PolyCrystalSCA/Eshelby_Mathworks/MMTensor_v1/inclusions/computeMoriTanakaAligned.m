% function [C_effective, C_matrix] = computeMoriTanakaAligned(C_matrix,
%     nu_matrix, concentrations, C_inclusions, aspect_ratios)
%
% Estimate using Mori-Tanaka's method the elasticity tensor of a
% heterogeneous medium with an isotropic matrix and axis-aligned
% inclusions.
%
% C_matrix: Young's modulus or stiffness tensor of the (isotropic) matrix
% nu_matrix: Poisson's ratio of the (isotropic) matrix
% concentrations: array with the concentrations of the inclusions
% C_inclusions: 6x6xN  array with the stiffness tensors of the inclusions
% aspect_ratios: either a Nx1 or a Nx3 array of the aspect ratios of the
% inclusions

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

function [C_effective, C_matrix] = computeMoriTanakaAligned(C_matrix, nu_matrix, concentrations, C_inclusions, aspect_ratios)

spheroidal = size(aspect_ratios,2) == 1;

% Compute matrix stiffness tensor:
if (prod(size(C_matrix)) == 1)
    % The matrix is assumed to be isotropic!
    C_matrix = makeIsotropicStiffnessTensor(C_matrix, nu_matrix);
end

% Compute matrix compliance tensor:
S_matrix = invertS4(C_matrix);

% Compute identity tensor
I = makeS4(1);

nb_inclusions = length(concentrations);
matrix_concentration = 1.0 - sum(concentrations);

% Compute the polarization(P) and strain concentration(T) tensors
% The (concentration weighted) polarisation(P) tensor
P_total = makeS4(0);

% The strain concentration(T) tensor
T_total = matrix_concentration*I;

for i = 1:nb_inclusions
  if (spheroidal)
    T = invertS4(I + makeSpheroidalEshelby(aspect_ratios(i), nu_matrix)*S_matrix*(C_inclusions(:,:,i)-C_matrix));
  else
    T = invertS4(I + makeEshelby(aspect_ratios(i,:), nu_matrix)*S_matrix*(C_inclusions(:,:,i)-C_matrix));
  end
  T_total = T_total + concentrations(i)*T;
  P_total = P_total + concentrations(i)*(C_inclusions(:,:,i)-C_matrix)*T;
end

% Compute the effective stiffness tensor
C_effective = C_matrix + P_total*invertS4(T_total);

end