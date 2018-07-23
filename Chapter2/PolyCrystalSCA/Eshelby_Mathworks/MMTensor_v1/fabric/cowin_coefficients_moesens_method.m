% function abc = cowin_coefficients_moesens_Method(E, A)
% 
% Compute Cowin's coefficients a, b, and c as defined in his 1985 paper. 
% This method uses a LU-based solution procedure proposed by Moesen.
% 
% E is a 6x6 symmetric matrix representing an orthogonal stiffness tensor
%   in sqrt2 notation. The coefficients of E should correspond to 
%   the principal axis system spanned by the eigenvectors of the fabric tensor. 
% A is a 3x1 vector containing the eigenvalues of the fabric tensor.
% 
% Return: abc is a 9x1 vector containing coefficients [a; b; c].

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

function abc = cowin_coefficients_moesens_Method(E, A)

%% Initialize fabric tensor A

% Make sure A is a 3x1 vector
A = reshape(A, 3, 1);

% Precompute the square of A
A2 = A.^2;

%% Prepare the RHS
AB_rhs = zeros(6,1);
AB_rhs(1:3) = diag(E(1:3,1:3));
AB_rhs(4:6) = [E(2,3); E(1,3); E(1,2)];

C_rhs = 0.5*diag(E(4:6,4:6));

%% Solve for the shear terms
C = [ 1   A(2)+A(3)   A2(2)+A2(3);
      1   A(1)+A(3)   A2(1)+A2(3);
      1   A(1)+A(2)   A2(1)+A2(2) ];

c = C\C_rhs;

%% Solve for the other terms 
AB(1:6, 1) = 1;
AB(1:3, 2) = 2*A;
AB(4:6, 2) = [A(2)+A(3); A(1)+A(3); A(1)+A(2)];
AB(1:3, 3) = 2*A2;
AB(4:6, 3) = [A2(2)+A2(3); A2(1)+A2(3); A2(1)+A2(2)];
AB(1:3, 4) = A2;
AB(4:6, 4) = [A(2).*A(3); A(1).*A(3); A(1).*A(2)];
AB(1:6, 5) = AB(:,2).*AB(:,4);
AB(1:6, 6) = AB(:,4).^2;

%%% Full solution used for debugging
%ABC = zeros(9,9);
%ABC(1:6,1:6) = AB;
%ABC(1:3,7:9) = 2*AB(1:3,1:3);
%ABC(7:9,7:9) = C;
%abc_2 = ABC\[AB_rhs; C_rhs];
  
AB_rhs(1:3) = AB_rhs(1:3) - 2*AB(1:3,1:3)*c;

ab = AB\AB_rhs;

%% Put everything together 
abc = [ab; c];

end