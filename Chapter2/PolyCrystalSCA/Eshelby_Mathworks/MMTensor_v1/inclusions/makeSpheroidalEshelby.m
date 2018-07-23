% function S = makeSpheroidalEshelby(alpha, nu)
%
% Compute the Eshelby tensor for axisymmetrical ellipsiodal inclusions with 
% aspect ratio alpha.
%
% Source: Mura(1987) p74-84
% nu: Poisson's ratio of the (isotropic) matrix material.
% alpha: the aspect ratio of the inclusions and has to be positive. Distinct cases are
%   penny-shape (alpha < 1), 
%   spheres (alpha = 1), 
%   cylinders (alpha = infty), and
%   oblate spheroid (alpha > 1).  
% A numerical instability of the implemented formulas must be remarked for alpha close but not equal to one.  

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

function S = makeSpheroidalEshelby(alpha, nu)

S = makeS4(0);

if (alpha < 0)
    return;
elseif (alpha == 1.0) % Spheres
    fnu = 1/(15*(1-nu));
    e_1111 = fnu*(7-5*nu);
    e_3333 = e_1111;
    e_1122 = fnu*(5*nu-1);
    e_1133 = e_1122;
    e_3311 = e_1122;
    e_1212 = fnu*(4-5*nu);
    e_1313 = fnu*(4-5*nu);
elseif (isinf(alpha)) % Needles
    hnu = 0.5/(1-nu);
    e_1111 = 0.25*hnu*(5-4*nu);
    e_3333 = 0;
    e_1122 = 0.25*hnu*(4*nu-1);
    e_1133 = nu*hnu;
    e_3311 = 0;
    e_1212 = 0.25*(3-4*nu)*hnu;
    e_1313 = 0.25;
else
  as = alpha^2;
  asm1 = as-1.0;
  hnu = 1/(2*(1-nu));

if (alpha < 1)
    q = (alpha/((-asm1)^(3/2)))*(acos(alpha) - alpha*sqrt(-asm1));
else
    q = (alpha/(asm1^(3/2)))*(alpha*sqrt(asm1) - acosh(alpha));
end

e_1111 = 0.75*hnu*(as/asm1)+ 0.5*hnu*(1-2*nu-9/4/asm1)*q;
e_3333 = hnu*(1-2*nu+(3*as-1)/asm1-(1-2*nu + 3*as/asm1)*q);
e_1122 = 0.5*hnu*(0.5*as/asm1-(1-2*nu + 0.75/asm1)*q);
e_1133 = hnu*(-as/asm1+0.5*(3*as/asm1 - (1-2*nu))*q);
e_3311 = hnu*(2*nu-1-1/asm1+(1-2*nu + 1.5/asm1)*q);
e_1212 = 0.5*hnu*(0.5*as/asm1 + (1-2*nu - 0.75/asm1)*q);
e_1313 = 0.5*hnu*(1-2*nu-(as+1)/asm1 - 0.5*(1-2*nu - 3*(as+1)/asm1)*q);
end

S = setS4Element(S,1,1,1,1, e_1111);
S = setS4Element(S,2,2,2,2, e_1111);
S = setS4Element(S,3,3,3,3, e_3333);

S = setS4Element(S,1,1,2,2, e_1122);
S = setS4Element(S,2,2,1,1, e_1122);

S = setS4Element(S,1,1,3,3, e_1133);
S = setS4Element(S,2,2,3,3, e_1133);

S = setS4Element(S,3,3,1,1, e_3311);
S = setS4Element(S,3,3,2,2, e_3311);

S = setS4Element(S,1,2,1,2, e_1212);
S = setS4Element(S,1,3,1,3, e_1313);
S = setS4Element(S,2,3,2,3, e_1313);
