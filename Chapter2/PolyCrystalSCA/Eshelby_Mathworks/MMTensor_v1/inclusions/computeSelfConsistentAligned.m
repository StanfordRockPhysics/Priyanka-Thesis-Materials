% function C_effective = computeSelfConsistentAligned(C_matrix, nu_matrix, concentrations, C_inclusions, aspect_ratio, max_nb_iterations, tolerance)
%
% Estimate using Hill's self consistent method the elasticity tensor of a
% heterogeneous medium with an isotropic matrix and axis-aligned
% inclusions.
%
% C_matrix: Young's modulus or stiffness tensor of the (isotropic) matrix
% nu_matrix: Poisson's ratio of the (isotropic) matrix
% concentrations: array with the concentrations of the inclusions
% C_inclusions: 6x6xN array with the stiffness tensors of the inclusions
% aspect_ratio: either a Nx1 or a Nx3 array of the aspect ratios of the
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


function C_effective = computeSelfConsistentAligned(C_matrix, nu_matrix, concentrations, C_inclusions, aspect_ratios, max_nb_iterations, tolerance)

spheroidal = size(aspect_ratios,2) == 1;

% Compute matrix stiffness tensor:
if numel(C_matrix) == 1
    % The matrix is assumed to be isotropic!
    C_matrix = makeIsotropicStiffnessTensor(C_matrix, nu_matrix);
end

% Initialize the effective stiffness with the matrix stiffness
C_effective = C_matrix;

% Add the matrix to the inclusions
concentrations = [ (1.0-sum(concentrations)) concentrations ];
C_inclusions = cat(3,C_matrix,C_inclusions);
aspect_ratios = [ones(1,size(aspect_ratios, 2)); aspect_ratios];

% Compute identity tensor
I = makeS4(1);

% The number of phases (matrix + inclusions)
nb_inclusions = length(concentrations);

% Initialize the polarization(P) and strain concentration(T) tensors
P_total = makeS4(0);
T_total = makeS4(0);

% This is actually not correct as the eshelby tensor should be computed
% based on the effective medium and not the matrix material
S_eshelby = zeros(6,6,nb_inclusions);
if (spheroidal)
  for i = 1:nb_inclusions
    S_eshelby(:,:,i) = makeSpheroidalEshelby(aspect_ratios(i), nu_matrix);
  end
else
  for i = 1:nb_inclusions
    S_eshelby(:,:,i) = makeEshelby(aspect_ratios(i,:), nu_matrix);
  end
end

for i = 1:max_nb_iterations
  % Compute matrix compliance tensor:
  S_effective = invertS4(C_effective);

for j = 1:nb_inclusions
  T = invertS4(I + S_eshelby(:,:,j)*S_effective*(C_inclusions(:,:,j)-C_effective)); 
  %T = invertS4(I + makeSpheroidalEshelby(aspect_ratio, nu_matrix)*S_effective*(C_inclusions{j}-C_effective));
  T_total = T_total + concentrations(j)*T;
  P_total = P_total + concentrations(j)*C_inclusions(:,:,j)*T;
end

  % Compute the effective stiffness tensor
  C_old = C_effective;
  C_effective = P_total*invertS4(T_total);
  abs_error = max(max(abs(C_effective-C_old)));
  
  if (abs_error < tolerance)
      fprintf('Convergence: absolute error after %d iterations: %g\n', i, abs_error);
      return
  end

end

fprintf('Maximum number of iterations reached!\nAbsolute error after %d iterations: %g\n', max_nb_iterations, abs_error);