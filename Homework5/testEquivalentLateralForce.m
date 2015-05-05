% Chopra example
% Miranda MDOF Seismic Analysis A, Page 14
Tn = 0.74 %s
Ie = 1.0;
Cd = 5.5;
Hi = 12; %ft
Cs = 0.105; %g
mass = 100; % kips/g
stiffness = 31.54; %kips/in
nfloors = 5;
[M,K] = computeMatrices(nfloors,mass,stiffness);
equivalentLateralForce(Tn,Cs,Cd,Ie,M,K,Hi)