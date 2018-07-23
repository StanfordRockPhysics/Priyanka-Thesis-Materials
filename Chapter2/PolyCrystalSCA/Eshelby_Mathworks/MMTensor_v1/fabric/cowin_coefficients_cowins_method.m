% function abc = cowin_coefficients_cowins_method(E, A)
% 
% Compute Cowin's coefficients a, b, and c as defined in his 1985 paper. 
% This method uses a quasi-direct solution procedure proposed by himself.
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


function abc = cowin_coefficients_cowins_method(E, A)

%% Initialize fabric tensor A

% Keep only the diagonal, the other coefficients should be zero.
if isequal(size(A),[3,3])
  A = diag(A);
end

% Make sure A is a 3x1 vector
A = reshape(A, 3, 1);

% Precompute the square of A
A2 = A.^2;

%% Initialize stiffness tensor E
E_diag = diag(E(1:3,1:3));
E_off = [E(2,3); E(1,3); E(1,2)]; % Off-diagonal terms
E_shear = 0.5*diag(E(4:6,4:6)); % Division by 2 due to sqrt2 notation. 

%% Initialize the coefficient matrices 
Z = [ ones(3,1) 2*A  2*A2];
B = repmat(A2, 1, 3) .* [ ones(3,1) 2*A  A2];
C = [ 1   A(2)+A(3)   A2(2)+A2(3);
      1   A(1)+A(3)   A2(1)+A2(3);
      1   A(1)+A(2)   A2(1)+A2(2) ];
    
D = [ A(2)*A(3); A(1)*A(3); A(1)*A(2) ];
D = [ D C(:,2).*D D.*D ];

%% Solve the system according to Cowin's proposal  

C_inv = inv(C);

F = B - Z*C_inv*D;
F_inv = inv(F);

c = C_inv * E_shear;
b = F_inv*E_diag - F_inv*Z*C_inv*[E_off+2*E_shear];
a = C_inv * (E_off - D*b); 

%% Put everything together
abc = [a;b;c];

end