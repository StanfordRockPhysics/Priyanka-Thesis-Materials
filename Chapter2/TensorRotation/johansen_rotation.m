function [R] = johansen_rotation(theta,xi,phi)
%% function to compute rotation matrix 'R' when rotation angles theta, xi and phi are specified in radians.

%% theta, xi and phi defined as per Johansen et. al., 2004 (Effect of grain scale alignment on seismic anisotropy
% and reflectivity of shales)

%% Coded by Priyanka Dutta, Dec, 2014

%% computing the rotation matrix
A(1,1:3) = [cos(phi)*cos(xi)-sin(phi)*cos(theta)*sin(xi) cos(phi)*sin(xi)+sin(phi)*cos(theta)*cos(xi) sin(theta)*sin(phi)];
A(2,1:3) = [-sin(phi)*cos(xi)-cos(phi)*cos(theta)*sin(xi) -sin(phi)*sin(xi)+cos(phi)*cos(theta)*cos(xi) sin(theta)*cos(phi)];
A(3,1:3) = [sin(theta)*sin(xi) -sin(theta)*cos(xi) cos(theta)];
R = A.';
end