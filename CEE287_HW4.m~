%% CEE 287: Homework 4
% Max Ferguson
% Nicole Paul

%% PART A - Problem 1

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

%% Part A - Problem 2

% Structure information
Tn = 1.0;
E = 0.05;
alpha = -0.15;
u0 = 0;
v0 = 0;

% Loading GM Information
FileList = dir('LomaPrietaGM');
N = size(FileList, 1);
for i = 4:N
    load(fullfile('LomaPrietaGM',FileList(i).name));
end
convtog = 0.00101971621;
GM = convtog*[Parking(:,2) SLAC_1(:,2) SLAC_2(:,2) VA_1(:,2) VA_2(:,2)];
dt = Parking(2,1); % same for all GM
nGM = 5;

% Using root-finding algorithm to find Cy for mu=1 for each GM
fn = @(C, A) FindMu2(Tn, E, A, dt, u0, v0, C, alpha) - 1;
Ce = NaN(nGM, 1);
for i = 1:nGM
    fn_A = @(C) fn(C, GM(:,i));
    Ce(i) = fzero(fn_A, 0.5);
end

% Plotting entire curves
m = 100;
Cy_range = linspace(0.05,1,m);
mu_mat = zeros(m, nGM);
for i = 1:nGM
    for j = 1:m
        mu_mat(j,i) = FindMu2(Tn, E, GM(:,i), dt, u0, v0, Cy_range(j), alpha);
    end
end

% Finding when collapse occurs
horz_tol = 1.5;
horz_diff = -(mu_mat(2:end,:) - mu_mat(1:end-1,:));
Cc = zeros(nGM,1);
ind_c = zeros(nGM,1);
mu_c = zeros(nGM,1);
for i = 1:nGM
    inds_c = find((horz_diff(:,i))>horz_tol);
    ind_c(i) = inds_c(end);
    Cc(i) = Cy_range(ind_c(i));
    mu_c(i) = mu_mat(ind_c(i),i);
end

% Plotting
figure;
colorcell = {'r','m','g','b','c'};
for i = 1:nGM
    h(i) = plot(mu_mat(:,i),Cy_range,strcat(colorcell{i},'-')); hold on;
    plot(1,Ce(i),strcat(colorcell{i},'^'),'MarkerSize',8);
    plot(mu_c(i),Cc(i),strcat(colorcell{i},'v'),'MarkerSize',8);
end
grid on;
xlabel('Ductility Demand, \mu');
ylabel('Seismic Coefficient, C_y');
title('Bilinear SDOF: T_n=1.0s, \xi=5%, \alpha=-0.15 - Loma Prieta');
legend(h,'Parking','SLAC_1','SLAC_2','VA_1','VA_2','Location','best');
xlim([0 20]);
hold off;



