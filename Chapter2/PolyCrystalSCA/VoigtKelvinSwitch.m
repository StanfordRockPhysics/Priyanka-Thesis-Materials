function [ ctrans ] = VoigtKelvinSwitch( c, opt )

% VoigtKelvinSwitch: Switch from Viogt to Kelvin notation and vice versa.

% Inputs
% c: 6X6 matrix: Elastic stiffness/compliance tensor of single crystal in Voigt/Kelvin notation
% opt: scalar that can take values 1 or 2: 
% opt = 1 implies either input 'c' is Voigt notation stifness and will be transformed to Kelvin notation stiffness or, 'c' is in Kelvin notation compliance and will be transformed to Voigt notation compliance
% opt = 2 implies input 'c' is in Kelvin (Voigt) notation stiffeness (compliance) and will be transformed to Voigt (Kelvin) notation stiffness (compliance).

% Outputs

% Detailed explanation goes here

% Coded by Priyanka Dutta, 2017

%* Check for valid 'opt'
if ~ismember(opt,[1,2]) 
    error('Check input: opt should be a number with value 1 or 2 only');
end

%* Voigt to Kelvin
v2k = [1 1 1 sqrt(2) sqrt(2) sqrt(2);1 1 1 sqrt(2) sqrt(2) sqrt(2);1 1 1 sqrt(2) sqrt(2) sqrt(2);...
    sqrt(2) sqrt(2) sqrt(2) 2 2 2;sqrt(2) sqrt(2) sqrt(2) 2 2 2;sqrt(2) sqrt(2) sqrt(2) 2 2 2];
if opt == 1
    ctrans = v2k.*c;
end

%* Kelvin to Voigt
k2v = [1 1 1 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);1 1 1 1/sqrt(2) 1/sqrt(2) 1/sqrt(2);...
    1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1/2 1/2 1/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1/2 1/2 1/2;1/sqrt(2) 1/sqrt(2) 1/sqrt(2) 1/2 1/2 1/2];
if opt == 2
    ctrans = k2v.*c;
end

end

