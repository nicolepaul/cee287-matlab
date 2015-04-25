%% CEE 287: Homework 4
% Max Ferguson
% Nicole Paul

%% PART A

Tn = 1.5;
E = 0.05;
alpha = -0.15;
u0 = 0;
v0 = 0;
load(fullfile('LomaPrietaGM','Parking.mat'));
A = Parking(:,2)*0.00102;
dt = Parking(2,1);
Cy = 0.01;

[u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAvgAccAlpha_Cy(Tn, E, A, dt, u0, v0, Cy, alpha);