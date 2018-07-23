% function [A, B] = computeEshelbyIntegrals(a)
% 
% Numerically compute the integrals I1, I2, I3, I11, I12, ... I33 required 
% for computing Eshelby's tensor. The computation is based on Mura(1987), 
% But the integrals are not multiplied by a factor 2*pi
% 
% This function provides a straightforward (inefficient) implementation. 
% Its advantage is that it works for arbitrary strictly positive aspect
% ratios.
%
% a: a vector of length three containing the aspect ratios of the
% ellipsoid.
%
% A: 3x1 vector containing integrals I_i.
% B: 3x3 matrix containing integrals I_ij.

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

function [A, B] = computeEshelbyIntegrals(a)

FD = @(s) sqrt((a(1).^2+s).*(a(2).^2+s).*(a(3).^2+s));
F1 = @(s) 1./((a(1).^2+s).*FD(s));
F2 = @(s) 1./((a(2).^2+s).*FD(s));
%F3 = @(s) 1./((a(3).^2+s).*FD(s));

F11 = @(s) 1./((a(1).^2+s).*(a(1).^2+s).*FD(s));
F12 = @(s) 1./((a(1).^2+s).*(a(2).^2+s).*FD(s));
F13 = @(s) 1./((a(1).^2+s).*(a(3).^2+s).*FD(s));
%F21 = @(s) 1./((a(2).^2+s).*(a(1).^2+s).*FD(s));
F22 = @(s) 1./((a(2).^2+s).*(a(2).^2+s).*FD(s));
F23 = @(s) 1./((a(2).^2+s).*(a(3).^2+s).*FD(s));
%F31 = @(s) 1./((a(3).^2+s).*(a(1).^2+s).*FD(s));
%F32 = @(s) 1./((a(3).^2+s).*(a(2).^2+s).*FD(s));
F33 = @(s) 1./((a(3).^2+s).*(a(3).^2+s).*FD(s));

c = prod(a);

A = zeros(3,1);
B = zeros(3,3);

A(1)=  c*quadgk(F1, 0, Inf);
A(2) = c*quadgk(F2, 0, Inf);
A(3) = 2 - A(1) - A(2); %quadgk(F3, 0, Inf);

%B = c*[ quadgk(F11, 0, Inf) quadgk(F12, 0, Inf) quadgk(F13, 0, Inf); 
%      quadgk(F21, 0, Inf) quadgk(F22, 0, Inf) quadgk(F23, 0, Inf); 
%      quadgk(F31, 0, Inf) quadgk(F32, 0, Inf) quadgk(F33, 0, Inf)];

B(1,1) = c*quadgk(F11, 0, Inf);
B(2,2) = c*quadgk(F22, 0, Inf);
B(3,3) = c*quadgk(F33, 0, Inf);
B(1,2) = c*quadgk(F12, 0, Inf);
B(1,3) = c*quadgk(F13, 0, Inf);
B(2,3) = c*quadgk(F23, 0, Inf);
B(2,1) = B(1,2);
B(3,1) = B(1,3);
B(3,2) = B(2,3);

end