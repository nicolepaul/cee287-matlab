%% CEE 287: Homework 4
% Max Ferguson
% Nicole Paul

%% PART A

Tn = 1.5;
E = 0.05;
alpha = -0.15;
u0 = 0;
v0 = 0;
dt = 0.02;
dtnew = 0.01;
load(fullfile('OtherData','ElCentroGMTH.mat'));
A = A_ElCentro;
A = interpolateGM(A,dt,dtnew);
Cy = 0.129;

[u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAvgAccAlpha_Cy(Tn, E, A, dtnew, u0, v0, Cy, alpha);

%% 

load(fullfile('OtherData','ElCentro10s'));
Tn = 0.5;
E = 0;
dt = 0.002;
uy = 0.4182;
u0 = 0;
v0 = 0;
M = 1; % assumed, should cancel out
g = 386;
wn = 2*pi/Tn;
C = 2*E*wn*M;
K = wn^2*M;
Cy = uy*K/(M*g);
alpha = 0;
%[u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAvgPrev(Tn, E, A_ElCentro, dt, u0, v0, Cy, 0);
[u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAvgAccAlpha_Cy(Tn, E, A_ElCentro, dt, u0, v0, Cy, alpha);
