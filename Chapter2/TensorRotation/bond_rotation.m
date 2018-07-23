function [ OP ] = bond_rotation( IP, R, choice )

% bond_rotation: co-ordinate rotation of elastic stiffness/compliance
% tensors in Voigt notation

%% IP = input stiffness/compliance tensor in Voigt notation (6X6 matrix) in
% refernce co-ordinate system
%% R = rotation matrix such that op = R*ip where 'op' is a vector in the final
% co-ordinate system and 'ip' is a vector in the refence co-ordinate system
%% OP = output stiffness/compliance tensor in Voigt notation (6X6 matrix) in
% final co-ordinate system
%% choice = 1 for stiffness tensor, 2 for compliance

%% Reference: RP Handbook, 2nd Edition, Pg: 19-21 

%% Coded by Priyanka Dutta, Dec, 2014

%% Bond transform
if choice == 1
    M = [R(1,1)^2 R(1,2)^2 R(1,3)^2 2*R(1,2)*R(1,3) 2*R(1,3)*R(1,1) 2*R(1,1)*R(1,2);...
         R(2,1)^2 R(2,2)^2 R(2,3)^2 2*R(2,2)*R(2,3) 2*R(2,3)*R(2,1) 2*R(2,1)*R(2,2);...
         R(3,1)^2 R(3,2)^2 R(3,3)^2 2*R(3,2)*R(3,3) 2*R(3,3)*R(3,1) 2*R(3,1)*R(3,2);...
         R(2,1)*R(3,1) R(2,2)*R(3,2) R(2,3)*R(3,3) R(2,2)*R(3,3)+R(2,3)*R(3,2) R(2,1)*R(3,3)+R(2,3)*R(3,1) R(2,2)*R(3,1)+R(2,1)*R(3,2);...
         R(3,1)*R(1,1) R(3,2)*R(1,2) R(3,3)*R(1,3) R(1,2)*R(3,3)+R(1,3)*R(3,2) R(1,1)*R(3,3)+R(1,3)*R(3,1) R(1,1)*R(3,2)+R(1,2)*R(3,1);...
         R(1,1)*R(2,1) R(1,2)*R(2,2) R(1,3)*R(2,3) R(2,2)*R(1,3)+R(1,2)*R(2,3) R(1,1)*R(2,3)+R(1,3)*R(2,1) R(2,2)*R(1,1)+R(1,2)*R(2,1)];
    OP = M*IP*M.';
elseif choice == 2
    N = [R(1,1)^2 R(1,2)^2 R(1,3)^2 R(1,2)*R(1,3) R(1,3)*R(1,1) R(1,1)*R(1,2);...
         R(2,1)^2 R(2,2)^2 R(2,3)^2 R(2,2)*R(2,3) R(2,3)*R(2,1) R(2,1)*R(2,2);...
         R(3,1)^2 R(3,2)^2 R(3,3)^2 R(3,2)*R(3,3) R(3,3)*R(3,1) R(3,1)*R(3,2);...
         2*R(2,1)*R(3,1) 2*R(2,2)*R(3,2) 2*R(2,3)*R(3,3) R(2,2)*R(3,3)+R(2,3)*R(3,2) R(2,1)*R(3,3)+R(2,3)*R(3,1) R(2,2)*R(3,1)+R(2,1)*R(3,2);...
         2*R(3,1)*R(1,1) 2*R(3,2)*R(1,2) 2*R(3,3)*R(1,3) R(1,2)*R(3,3)+R(1,3)*R(3,2) R(1,1)*R(3,3)+R(1,3)*R(3,1) R(1,1)*R(3,2)+R(1,2)*R(3,1);...
         2*R(1,1)*R(2,1) 2*R(1,2)*R(2,2) 2*R(1,3)*R(2,3) R(2,2)*R(1,3)+R(1,2)*R(2,3) R(1,1)*R(2,3)+R(1,3)*R(2,1) R(2,2)*R(1,1)+R(1,2)*R(2,1)];
    OP = N*IP*N.';
else
    error('Check input: choice');
end

