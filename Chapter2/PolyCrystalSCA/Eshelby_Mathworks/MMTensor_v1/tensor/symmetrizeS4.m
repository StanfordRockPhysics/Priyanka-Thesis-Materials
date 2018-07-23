% function S = symmetrizeS4(T)
%
% Compute the symmetrization of symmetric fourth order tensor T.
% This symmetrization goes much further than just computing 0.5*(T+T').
% It symmetrizes the elements such that 
% T_{ijkl} = sum_{ijkl} / 24 
% {ijkl} are all 24 = 4! permutations {ijkl, ijlk, ikjl, iklj, ...}

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

function S = symmetrizeS4(T)
S = makeS4(0);
SQRT2 = sqrt(2);

%% Compute the symmetrization of the following part

%  A  F  E  0  0  0 
%  F  B  D  0  0  0
%  E  D  C  0  0  0
%  0  0  0  D  0  0
%  0  0  0  0  E  0
%  0  0  0  0  0  F

% Terms T(1,1) T(2,2) and T(3,3) remain the same
S(1,1) = T(1,1);
S(2,2) = T(2,2);
S(3,3) = T(3,3);

% Terms T(2,1) T(1,2) and T(6,6)
S(2,1) = (T(2,1)+T(1,2)+2*T(6,6)) / 6;
S(1,2) = S(2,1);
S(6,6) = 2*S(2,1);

% Terms T(3,1) T(1,3) and T(5,5)
S(3,1) = (T(3,1)+T(1,3)+2*T(5,5)) / 6;
S(1,3) = S(3,1);
S(5,5) = 2*S(3,1);

% Terms T(3,2) T(2,3) and T(4,4)
S(3,2) = (T(3,2)+T(2,3)+2*T(4,4)) / 6;
S(2,3) = S(3,2);
S(4,4) = 2*S(3,2);

%% Now compute the symmetrization of 
%  *  *  *  A  0  0 
%  *  *  *  0  B  0
%  *  *  *  0  0  C
%  A  0  0  *  C  B
%  0  B  0  C  *  A
%  0  0  C  B  A  *
% The asterisk terms (*) were already done before

% Term A
S(4,1) = (T(4,1)+T(1,4)+SQRT2*(T(6,5)+T(5,6))) / 6;
S(1,4) = S(4,1);
S(5,6) = SQRT2*S(4,1);
S(6,5) = S(5,6);

% Term B 
S(5,2) = (T(5,2)+T(2,5)+SQRT2*(T(6,4)+T(4,6))) / 6;
S(2,5) = S(5,2);
S(4,6) = SQRT2*S(5,2);
S(6,4) = S(4,6);

% Term C
S(6,3) = (T(6,3)+T(3,6)+SQRT2*(T(5,4)+T(4,5))) / 6;
S(3,6) = S(6,3);
S(4,5) = SQRT2*S(6,3);
S(5,4) = S(4,5);

%% Finally compute the symmetrization of the remaining terms
% These terms won't change if T is already symmetric such that T=T', 
% or alternatively T_ijkl = T_klij.
%  *  *  *  *  A  B 
%  *  *  *  F  *  C
%  *  *  *  E  D  *
%  *  F  E  *  *  *
%  A  *  D  *  *  *
%  B  C  *  *  *  *
% The asterisk terms (*) were already done before

S(5,1) = 0.5*(T(5,1)+T(1,5));
S(1,5) = S(5,1);

S(6,1) = 0.5*(T(6,1)+T(1,6));
S(1,6) = S(6,1);

S(6,2) = 0.5*(T(6,2)+T(2,6));
S(2,6) = S(6,2);

S(5,3) = 0.5*(T(5,3)+T(3,5));
S(3,5) = S(5,3);

S(4,3) = 0.5*(T(4,3)+T(3,4));
S(3,4) = S(4,3);

S(4,2) = 0.5*(T(4,2)+T(2,4));
S(2,4) = S(4,2);

%% Done!
end
