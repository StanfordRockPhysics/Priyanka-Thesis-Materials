% function S = makeEshelby(a, nu)
%
% Compute the Eshelby tensor for general ellipsoidal inclusions.
% Source: Mura(1987) p74-84 
%
% a: a vector of length three containing the aspect ratios of the
%    ellipsoidal inclusions.
% nu: Poisson's ratio of the (isotropic) matrix material.

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

function S = makeEshelby(a, nu)

[A,B] = computeEshelbyIntegrals(a);
a = a.^2;

c_offd = 1./(4*(1-nu)); % Off-diagonal coefficient
c_diag = 3*c_offd; % Diagonal coefficients
c_maj = (1-2*nu)*c_offd;

S = zeros(6,6);

S(1,1) = c_maj*A(1) + c_diag*a(1)*B(1,1);
S(2,2) = c_maj*A(2) + c_diag*a(2)*B(2,2);
S(3,3) = c_maj*A(3) + c_diag*a(3)*B(3,3);
S(1,2) = -c_maj*A(1) + c_offd*a(2)*B(1,2);  
S(1,3) = -c_maj*A(1) + c_offd*a(3)*B(1,3);
S(2,1) = -c_maj*A(2) + c_offd*a(1)*B(2,1);
S(2,3) = -c_maj*A(2) + c_offd*a(3)*B(2,3);
S(3,1) = -c_maj*A(3) + c_offd*a(1)*B(3,1);
S(3,2) = -c_maj*A(3) + c_offd*a(2)*B(3,2);
S(4,4) = (a(2)+a(3))*B(2,3)*c_offd + c_maj*(A(2)+A(3));
S(5,5) = (a(1)+a(3))*B(1,3)*c_offd + c_maj*(A(1)+A(3));
S(6,6) = (a(1)+a(2))*B(1,2)*c_offd + c_maj*(A(1)+A(2));

end