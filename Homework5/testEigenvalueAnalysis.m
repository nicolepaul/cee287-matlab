
% Test eigenvalueAnalysis based on the example in
% Miranda, MDOF Structures Part A, page 32
mass = 100/386.4; % At each floor
stiffness = 31.54; % At each floor
nfloors = 5;
[w,T,sphi,Gamma] = eigenvalueAnalysis(nfloors,mass,stiffness);

w
T
Gamma
printModalVectorTable(sphi);


