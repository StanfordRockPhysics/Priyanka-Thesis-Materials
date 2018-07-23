function plotDirectionalRoseDiagrams(P)


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

%[N2, N4] = fabric_moment_tensor(D)
%[D2,D4] = fabric_decomposition(N2, N4)
%[F2,F4] = fabric_tensor(N2,N4)
%[p2, p4] = fabric_fitness_test(nb_edges,D2,D4)

%% Plot the data points as data

[N2,N4] = fabric_moment_tensor(P)
[F2,F4] = fabric_tensor(N2,N4)

%% Plot the representation surface of the fabric tensor
[x,y,z] = sphere(50);
S = size(x);

x = reshape(x, 1, numel(x));
y = reshape(y, 1, numel(y));
z = reshape(z, 1, numel(z));

N = makeS2FromVector([x;y;z]);

sN2 = N2'*N;
sF2 = F2'*N;

sN4 = sN2;
sF4 = sF2;

for i = 1:length(sN4)
  sN4(i) = N(:,i)'*N4*N(:,i);
  sF4(i) = N(:,i)'*F4*N(:,i);
end

x = reshape(x, S);
y = reshape(y, S);
z = reshape(z, S);
sN2 = reshape(sN2, S);
sF2 = reshape(sF2, S);
sN4 = reshape(sN4, S);
sF4 = reshape(sF4, S);

size(P)
P_norm = euclidean_norm(P);
P_total = mean(P_norm);
P = [-P;P]/P_total;

%% Draw the star plot of the points
figure;
plot3(P([1 4],:), P([2 5],:), P([3 6],:), 'k.-');
colormap gray;
axis equal
axis vis3d
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('Star plot')
set(gca,'Color','none');
set(gcf,'Color','none');

%% Plot the representation surface of N2
figure;
surf(sN2.*x,sN2.*y,sN2.*z, sN2, 'FaceColor','none','EdgeColor','blue');%,'EdgeColor','none');%surf(X,Y,Z,'FaceColor','gray','EdgeColor','none')
colormap gray;
axis vis3d
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('Moment tensor N2')
set(gca,'Color','none');
set(gcf,'Color','none');

%% Plot the representation surface of N4
figure;
surf(sN4.*x,sN4.*y,sN4.*z, sN4, 'FaceColor','none','EdgeColor','red');%,'EdgeColor','none');%surf(X,Y,Z,'FaceColor','gray','EdgeColor','none')
colormap gray;
axis equal
axis vis3d
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('Moment tensor N4')
set(gca,'Color','none');
set(gcf,'Color','none');

%% Plot the representation surface of F2 and F4
figure;
plot3(P([1 4],:), P([2 5],:), P([3 6],:), 'k.'); hold on;
surf(sF2.*x,sF2.*y,sF2.*z, sF2, 'FaceColor','none','EdgeColor','blue');%,'EdgeColor','none');%surf(X,Y,Z,'FaceColor','gray','EdgeColor','none')
surf(sF4.*x,sF4.*y,sF4.*z, sF4, 'FaceColor','none','EdgeColor','red');%,'EdgeColor','none');%surf(X,Y,Z,'FaceColor','gray','EdgeColor','none')

ca = caxis;
ca = [min(0,ca(1)),ca(2)];
caxis(ca);
colormap gray;
axis equal
axis vis3d
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('Fabric tensors F2(blue) and F4(red)')
set(gca,'Color','none');
set(gcf,'Color','none');

end

