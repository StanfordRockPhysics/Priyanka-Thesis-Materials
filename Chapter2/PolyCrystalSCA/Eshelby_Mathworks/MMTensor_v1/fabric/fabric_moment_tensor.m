% function [N2, N4, w] = fabric_moment_tensor(p, w)
%
% Compute the fabric moment tensors N2 and N4 for the direction vectors
% in p, weighted by the weights in w.
%
% p: a 3 x N matrix with N direction vectors.
% The norm of the direction vectors is used for weighting.
% w: a 1 x N vector with N (optional) additional weighting factors.
% Returned w: weight vector corrected for the norms of the direction vectors.

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


function [N2, N4, w] = fabric_moment_tensor(p, w)

%% Define global constants and check the input arguments

SQRT2 = sqrt(2);

%if nargin == 0
%  n = make_uniform_directions(8)
%end 

if nargin <= 1
  w = ones(1,size(p,2));
end

%% Normalize the direction vectors
p_norm = euclidean_norm(p);
p = p ./ repmat(p_norm, 3, 1);

%% Normalize the weights
w = w .* p_norm;
w = w ./ sum(w);

%% Compute the moments in sqrt2 notation

% The inner products of the directions with themselves in sqrt2 notation 
NI = makeS2FromVector(p);

% Compute the weighted second order moment.
N2 = NI*w';

% Compute the weighted fourth order moment.
N4 = repmat(w,6,1).*NI*NI';

end 