% function oc = orthogonal_coefficients(C,Q,N)
%
% Compute the coefficients of the stiffness tensors in C with
% respect to the orthogonal basis Q.  
%
% C: 6 x 6 x n matrix containing n stiffness tensors in matrix notation.
%    In this case the non-orthogonal coefficients are discarded.
% C: 9 x n matrix with in each column a stiffness tensor in vector notation.  
% Q: 9 x r coefficient matrix of the rank r orthogonal basis. 
% N: 9 x r vector containing the squared norms of the orthogonal basis.
% The size r of the input vectors/matrices may be 2, 5 or 9 and can be
% larger than the actual rank of the orthogonal basis.
% 
% oc: r x n matrix with the orthogonal coefficients and 
%     r is the actual rank of the orthogonal basis 
% See also: orthogonal_basis_numerical, orthogonal_basis_analytical

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

function oc = orthogonal_coefficients(C,Q,N)

%% Convert the stiffness tensors to vector notation if necessary.
if size(C,1) == 6
  C = matrix_to_vector_notation(C);
end

%% Determine the rank and the sizes of the matrices
n = size(C,2); % Number of stiffness tensors
r = rank_from_norms(N); % Rank of the orthogonal basis

%% Compute the coefficients
oc = Q(:,1:r)'*C./repmat(N(1:r), [1 n]);

end


