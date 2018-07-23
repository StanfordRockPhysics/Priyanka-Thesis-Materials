% function plotRepresentationSurface(C)
%
% Plot the representation surface of a stiffness tensor C.
%
% C: 6x6 matrix representing a stiffness tensor in sqrt2 matrix notation.

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

function plotRepresentationSurface(C)

%% Assume that C has sqrt2 notation.

S = inv(C);

[x,y,z] = sphere(50);
c = zeros(size(z));

for i = 1:numel(x)
  %E = orthorhombic_modulus(S, [x(i); y(i); z(i)]);
  E = triclinic_modulus(S, [x(i); y(i); z(i)]);
  x(i) = E*x(i);
  y(i) = E*y(i);
  z(i) = E*z(i);
  c(i) = E;
end

%c = max(max(c)) - c;

figure;
surf(x,y,z, c);%,'EdgeColor','none');%surf(X,Y,Z,'FaceColor','gray','EdgeColor','none')
ca = caxis;
ca = [min(0,ca(1)),ca(2)];
caxis(ca);
colormap gray;
axis equal
axis vis3d
xlabel('x1');
ylabel('x2');
zlabel('x3');
set(gca,'Color','w');
set(gcf,'Color','w');

end


function E = orthorhombic_modulus(S, n)
n = n.^2;
N = n*n';

% The book of NYE uses the engineering stiffness tensor, not sqrt2 notation
% and S_eng(4:6,4:6) = 4*S_ijkl and S_sqrt2(4:6,4:6) = 2*S_ijkl.
% Thus S_eng = 2*S_sqrt(2), Hence the factor 2 for S(4,4), S(5,5), S(6,6).  
E = sum(sum(S(1:3, 1:3).*N)) + 2*(N(2,3)*S(4,4) + N(1,3)*S(5,5) + N(1,2)*S(6,6));
E = 1/E;
end

function E = triclinic_modulus(S, n)
N = makeS2FromVector(n);
% The book of NYE uses the engineering stiffness tensor, not sqrt2 notation
% and S_eng(4:6,4:6) = 4*S_ijkl and S_sqrt2(4:6,4:6) = 2*S_ijkl.
% Thus S_eng = 2*S_sqrt(2), Hence the factor 2 for S(4,4), S(5,5), S(6,6).
E = N'*S*N;

%Alternative computation: 
%N = N'*N;
%E = sum(sum(N.*S));

E = 1/E;
end

