% function d = eval_fabric_distribution_function(F, theta, phi)
%
% Integrate given fabric tensor F the distribution probability in the
% angular intervals given by theta and phi.
%
% The fabric tensor of the second kind in sqrt2 notation
% theta: The angle interval with the vertical axis, typically between [0, pi)
% phi: The angle interval in the horizontal plane with the x-axis in [0, 2*pi)

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

function d = eval_fabric_distribution_function(F, theta, phi)

fun = @(a,b) 1;

%% If F is a second order tensor
if size(F, 2) == 1
  fun = @(a, b) eval_integral_fun_S2(F,a,b);
elseif size(F,2) == 6
  fun = @(a, b) eval_integral_fun_S4(F,a,b);
end

d = dblquad(fun, theta(1), theta(2), phi(1), phi(2));
d = d/(4*pi);

end

function i = eval_integral_fun_S2(F, a, b)
  
sa = sin(a);
ca = cos(a);
sb = sin(b);
cb = cos(b);

n = [sa.*cb; sa.*sb; ones(size(a)).*ca];
fun_value = contractS2ToScalar(F, n);

i = fun_value.*sa;

end

function i = eval_integral_fun_S4(F, a, b)
  
sa = sin(a);
ca = cos(a);
sb = sin(b);
cb = cos(b);

n = [sa.*cb; sa.*sb; ones(size(a)).*ca];
N = makeS2FromVector(n);

fun_value = zeros(size(a));
for ai = 1:numel(a)
  fun_value(ai) = N(:,ai)'*F*N(:,ai);
end

i = fun_value.*sa;

end