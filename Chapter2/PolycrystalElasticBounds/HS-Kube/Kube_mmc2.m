%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to test bounds.m. Each case corresponds to a different material.
% Output C corresponds with the values given in Table 1. 
% All elastic constants given in GPa. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all
N = 50;

for i = 1:19
  i
    switch i
        case 1,
            % Elastic constants of Dibenzyl
            c11 = 9.81; c12 = 4.48; c13 = 4.25; c14 = 0; c15 = 2.91; c16 = 0; c22 = 7.37; c23 = 4.78; c24 = 0; c25 = 0.92; c26 = 0; c33 = 8.06; c34 = 0; c35 = 0.27; c36 = 0; c44 = 3.24; c45 = 0; c46 = -0.81; c55 = 2.72; c56 = 0; c66 = 2.62;
        case 2,
            % Elastic constants of Ethylene diamine tartrate
            c11 = 64.0; c12 = 13.0; c13 = 26.7; c14 = 0; c15 = 13.1; c16 = 0; c22 = 33.6; c23 = 10.3; c24 = 0; c25 = -0.1; c26 = 0; c33 = 22.5; c34 = 0; c35 = 7.4; c36 = 0; c44 = 5.22; c45 = 0; c46 = -0.1; c55 = 11.6; c56 = 0; c66 = 5.23;
        case 3,
            % Elastic constants of Oxalic acid dihydrate
            c11 = 21.7; c12 = 9.46; c13 = 12.18; c14 = 0; c15 = -0.39; c16 = 0; c22 = 11.6; c23 = 6.14; c24 = 0; c25 = 0.96; c26 = 0; c33 = 40.1; c34 = 0; c35 = -7.0; c36 = 0; c44 = 2.08; c45 = 0; c46 = 0.36; c55 = 8.09; c56 = 0; c66 = 9.16;
        case 4,
            % Elastic constants of Potassium cobalt cyanide
            c11 = 21.2; c12 = 2.1; c13 = 13.0; c14 = 0; c15 = 0; c16 = 0; c22 = 28.1; c23 = 0; c24 = 0; c25 = -7.4; c26 = 0; c33 = 34.0; c34 = 0; c35 = 0.1; c36 = 0; c44 = 10.3; c45 = 0; c46 = 3.3; c55 = 10.8; c56 = 0; c66 = 9.7;
        case 5,
            % Elastic constants of Tin difluoride
            c11 = 47.9; c12 = 9.3; c13 = 5.3; c14 = 0; c15 = -5.1; c16 = 0; c22 = 20.6; c23 = 14.8; c24 = 0; c25 = 3.1; c26 = 0; c33 = 33.6; c34 = 0; c35 = 6.5; c36 = 0; c44 = 17.4; c45 = 0; c46 = 3.5; c55 = 12.9; c56 = 0; c66 = 14.4;
        case 6,
            % Elastic constants of Ammonium tetroxalate dihydrate
            c11 = 21.9; c12 = 12.0; c13 = 10.4; c14 = 1.6; c15 = 6.0; c16 = -1.0; c22 = 45.9; c23 = 16.3; c24 = 11.6; c25 = 2.0; c26 = -3.8; c33 = 36.4; c34 = 3.8; c35 = 2.0; c36 = -0.8; c44 = 10.4; c45 = 0.1; c46 = 0.1; c55 = 5.40; c56 = 0.1; c66 = 4.4;
        case 7,
            % Elastic constants of Cesium trihydrogen selenite
            c11 = 57.1; c12 = 20.6; c13 = 31.6; c14 = -4.3; c15 = -0.4; c16 = -2.2; c22 = 35.8; c23 = 23.4; c24 = -2.8; c25 = -0.1; c26 = -0.6; c33 = 58.4; c34 = -0.8; c35 = -2.8; c36 = -0.8; c44 = 16.5; c45 = -1.8; c46 = 1.2; c55 = 15.2; c56 = -3.5; c66 = 12.0;
        case 8,
            % Elastic constants of Copper sulfate pentahydrate
            c11 = 56.5; c12 = 26.5; c13 = 32.1; c14 = -3.3; c15 = -0.8; c16 = -3.9; c22 = 43.3; c23 = 34.7; c24 = -0.7; c25 = -2.1; c26 = 0.2; c33 = 56.9; c34 = -4.4; c35 = -2.1; c36 = -1.6; c44 = 17.3; c45 = 0.9; c46 = 0.3; c55 = 12.2; c56 = -2.6; c66 = 10.0;
        case 9,
            % Elastic constants of Lithium hydrogen oxalate monohydrate
            c11 = 95.2; c12 = 26.0; c13 = 28.3; c14 = 2.5; c15 = 30.7; c16 = 10.6; c22 = 86.7; c23 = 14.2; c24 = -9.3; c25 = 4.7; c26 = 1.2; c33 = 35.0; c34 = -0.5; c35 = 9.2; c36 = 1.3; c44 = 9.6; c45 = 1.6; c46 = 4.3; c55 = 19.4; c56 = 0.9; c66 = 22.2;
        case 10,
            % Elastic constants of Potassium tetroxalate dihydrate
            c11 = 25.4; c12 = 11.8; c13 = 9.83; c14 = 0.72; c15 = 6.12; c16 = -1.23; c22 = 47.8; c23 = 14.0; c24 = 11.3; c25 = 1.46; c26 = -2.70; c33 = 34.3; c34 = 2.19; c35 = 1.47; c36 = 0.40; c44 = 10.2; c45 = -0.82; c46 = 0.53; c55 = 5.69; c56 = 0.70; c66 = 4.99;
        case 11,
            % Elastic constants of Sodium hydrogen oxalate hydrate
            c11 = 24.86; c12 = 16.81; c13 = 11.84; c14 = 1.36; c15 = -0.50; c16 = -0.60; c22 = 83.42; c23 = 20.85; c24 = -0.62; c25 = 3.30; c26 = 4.26; c33 = 78.78; c34 = -5.89; c35 = 1.42; c36 = 0.91; c44 = 15.37; c45 = 0.27; c46 = 2.46; c55 = 5.69; c56 = 1.03; c66 = 8.34;
        case 12,
            % Elastic constants of An0
            c11 = 68.335; c22 = 184.3432; c33 = 180.0076; c44 = 24.9772; c55 = 26.9008; c66 = 33.5526; c12 = 32.1813; c13 = 30.4224; c15 = -2.2533; c23= 4.9698; c25 = -7.7979; c35 = 7.485; c46 = -7.1946; c14 = 4.8676; c16 = -0.9297; c24 = -4.3754; c26 = -6.376; c34 = -9.1671; c36 = -9.3982; c45 = -2.4157; c56 = 0.6107;
        case 13,
            % Elastic constants of An25
            c11 = 87.0577; c22 = 174.9225; c33 = 166.0891; c44 = 22.909; c55 = 29.04; c66 = 35.0378; c12 = 43.9349; c13 = 35.3897; c15 = -0.401; c23= 17.9522; c25 = -2.8696; c35 = 4.5962; c46 = -5.1742; c14 = 6.132; c16 = -0.6003; c24 = -5.8733; c26 = -6.4679; c34 = -2.9144; c36 = -10.6912; c45 = -1.3141; c56 = 0.7947;
        case 14,
            % Elastic constants of An37
            c11 = 96.1748; c22 = 189.4431; c33 = 171.9478; c44 = 23.6214; c55 = 33.0783; c66 = 35.5199; c12 = 46.0683; c13 = 38.4058; c15 = -0.152; c23= 15.4201; c25 = -5.0887; c35 = 7.2247; c46 = -4.828; c14 = 5.8709; c16 = -0.4071; c24 = -6.9728; c26 = -6.7722; c34 = 2.2077; c36 = -9.8442; c45 = -1.1262; c56 = 1.3771;
        case 15,
            % Elastic constants of An48
            c11 = 104.6397; c22 = 201.3603; c33 = 172.7897; c44 = 22.9499; c55 = 33.0057; c66 = 35.6241; c12 = 51.5; c13 = 43.8532; c15 = 0.1319; c23= 14.4585; c25 = -4.8006; c35 = 6.9137; c46 = -3.8397; c14 = 6.4965; c16 = -0.801; c24 = -2.4317; c26 = -9.8523; c34 = -0.3603; c36 = -5.7297; c45 = -1.0003; c56 = 2.0575;
        case 16,
            % Elastic constants of An60
            c11 = 109.3177; c22 = 185.4687; c33 = 164.1266; c44 = 22.1677; c55 = 33.1172; c66 = 36.8025; c12 = 53.0594; c13 = 42.1414; c15 = 1.2076; c23= 21.9334; c25 = 0.7258; c35 = 2.4989; c46 = 1.3602; c14 = 7.6048; c16 = -7.7187; c24 = -2.9306; c26 = -6.8497; c34 = 0.244; c36 = 0.6965; c45 = 0.1938; c56 = 2.8167;
        case 17,
            % Elastic constants of An78
            c11 = 120.4077; c22 = 191.6068; c33 = 163.6588; c44 = 23.252; c55 = 32.7587; c66 = 35.0393; c12 = 56.5598; c13 = 49.8835; c15 = 3.2471; c23= 26.2717; c25 = 5.3709; c35 = 1.676; c46 = 0.927; c14 = 9.0466; c16 = -3.0004; c24 = 2.1181; c26 = -9.8744; c34 = 1.714; c36 = -8.0794; c45 = 0.7831; c56 = 4.5015;
        case 18,
            % Elastic constants of An96
            c11 = 132.2194; c22 = 200.2049; c33 = 163.9226; c44 = 24.6175; c55 = 36.5852; c66 = 35.9933; c12 = 63.9587; c13 = 55.2934; c15 = 5.0678; c23= 31.8699; c25 = 3.5499; c35 = 0.4826; c46 = -2.2449; c14 = 9.5249; c16 = -10.7509; c24 = 7.4737; c26 = -7.2259; c34 = 6.6292; c36 = 1.6277; c45 = 2.9544; c56 = 5.1855;
        case 19,
            % Elastic constants of Graphite
            c11 = 1060; c12 = 180; c13 = 15.0; c14 = 0; c15 = 0; c16 = 0; c22 = c11; c23 = c13; c24 = 0; c25 = 0; c26 = 0; c33 = 36.5; c34 = 0; c35 = 0; c36 = 0; c44 = 4.0; c45 = 0; c46 = 0; c55 = c44; c56 = 0; c66 = (1/2)*(c11-c12);
            
    end
    
    % Run bounds.m for each material
    eu = bounds(1,N,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66);
    el = bounds(2,N,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66);
    ou = bounds(3,N,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66);
    ol = bounds(4,N,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66);
    
    % Create matrix of the Voigt, Reuss, Self-consistent, and
    % Hashin-Shtrikman moduli
    C(i,:) = [ol(1,2) el(1,2) ou(1,N) eu(1,2) ou(1,2) ol(2,2) el(2,2) ou(2,N) eu(2,2) ou(2,2) ou(2,2)/ol(2,2)+5*ou(1,2)/ol(1,2)-6];
end
round(10*C)/10

