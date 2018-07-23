% function A = makeTensorTransformationMatrix(a)
%
% Construct from axis system a rotation tensor in matrix notation.
% 
% a: a 3 x 3 matrix of which the columns represent an axis system (i.e. a
%    rotation from new coordinates to old coordinates.)
% A: a 6 x 6 matrix representing a second order rotation matrix. 

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

function A = makeTensorTransformationMatrix(a)

s = sqrt(2);
p2 = [2 3 1];
p3 = [3 1 2];
A = zeros(6,6);
A(1:3,1:3) = a.*a;
A(1:3,4) = s*a(:,2) .* a(:,3);
A(1:3,5) = s*a(:,1) .* a(:,3);
A(1:3,6) = s*a(:,1) .* a(:,2);
A(4,1:3) = s*a(2,:) .* a(3,:);
A(5,1:3) = s*a(3,:) .* a(1,:);
A(6,1:3) = s*a(1,:) .* a(2,:);
A(4:6,4) = a(p2,2).*a(p3,3) + a(p2,3).*a(p3,2);
A(4:6,5) = a(p2,3).*a(p3,1) + a(p2,1).*a(p3,3);
A(4:6,6) = a(p2,1).*a(p3,2) + a(p2,2).*a(p3,1);