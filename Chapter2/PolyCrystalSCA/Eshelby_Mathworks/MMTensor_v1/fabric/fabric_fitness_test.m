% function [p_2, p_4] = fabric_fitness_test(N, D2, D4)
%
% Return the p-values for the non-uniformity of D2 and D4. Thus, p_2 and p_4
% are close to zero when D2 or D4 is uniform.
% N is the number of direction vectors in the directional data set. 

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

function [p_2, p_4] = fabric_fitness_test(N, D2, D4)

if N < 20
    fprintf('N (here: %d) should be large for fitness tests\n', N);
end

test_statistic_2 = 2*N/15*(D2'*D2);
p_2 = 1-chi2cdf(test_statistic_2,2);
fprintf('p_2 < 0.005 indicates a significant difference from a uniform distribution, so consider D2\n');
                           
if nargout < 2
    return
end

if nargin < 3
    p_2 = 1;
    return
end

D2norm = norm(D2);
if D2norm > 0.1
    fprintf('For testing second order fitness, D2(norm=%g) should be small in comparison to 1\n', D2norm);
end

% The test requires the tensor contraction Dijkl Dijkm Dlm
% Compute it the ugly way...
test_statistic_4 = 0;
D = D4*D4;
for k=1:3
    for l=1:3
        for m = 1:3
            test_statistic_4 = test_statistic_4 + getS4Element(D,k,l,k,m)*getS2Element(D2,l,m);
        end
    end
end

% Compute the test statistic
test_statistic_4 = 8*N/315*(sum(sum(D4.^2))-8/11*test_statistic_4);

% Compute the p-value
p_4 = 1-chi2cdf(test_statistic_4,9);

fprintf('p_2 < 0.005 indicates that fourth order terms should be considered.\n');

end
