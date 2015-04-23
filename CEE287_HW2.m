%% CEE 287: Homework 2
% Nicole Paul
% Max Ferguson

%% Part A
close all;

if (1)
    % Compare to Chopra pg 279
    load ElCentro10s
    Tn = 0.5;
    E = 0;
    dt = 0.002;
    uy = 0.4182;
    u0 = 0;
    v0 = 0;
    
    [u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAverageAcceleration(Tn, E, A_ElCentro, dt, u0, v0, uy);
    fprintf('Max displacement should be: 1.71 in\n')
    fprintf('Max displacement is: %.2f in\n',max(abs(u)))
    
else
    % Chopra Page 193
    % Must also change Newmark...Acceleration mass to 0.2533
    dt = 0.1;
    uy = 0.75;
    Tn = 1;
    E = 0.05;
    u0 = 0;
    v0 = 0;
    
    %M = 0.253; % assumed, should cancel out
    g = 386; % assumed input is in [g]
    GM = -[0, 5, 8.6603 10.000, 8.6603 5 0 0 0 0 0]'/(0.2533*g);
    NewmarkAverageAcceleration(Tn, E, GM, dt, u0, v0, uy);
end

%% Part B

load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
Tn = 0.4;
E = 0.05;
u0 = 0;
v0 = 0;
%%
fn = @(u) FindMu(Tn, E, A_Parking, dt, u0, v0, u) - 1;
M = 1; % assumed, should cancel out
g = 386;
wn = 2*pi/Tn;
C = 2*E*wn*M;
K = wn^2*M;

Cy = K*fzero(fn, 500/K)/(M*g);

Ce = 0.01:0.01:1;
uy = Ce*M*g/K;
n = numel(uy);
muvec = zeros(n,1);
for i = 1:n
    muvec(i) = FindMu(Tn, E, A_Parking, dt, u0, v0, uy(i));
end

figure;
plot(muvec,Ce,'b-', [1 1], [min(Ce) max(Ce)], 'r--', 1, Cy, 'mo'); grid on;
xlabel('\mu'); ylabel('C_y');
xlim([0 6]); ylim([min(Ce) max(Ce)]);
title('3-Story T_n = 0.4s - Loma Prieta - Parking Record');
legend('Calculated \mu','\mu = 1', 'Calculated root');

%%

load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
u0 = 0;
v0 = 0;
Cy = [0.6875 0.3429];
Tn = [0.4 1.0];
E = 0.05;

Sd = zeros(2,1);
Sa = zeros(2,1);
for i = 1:2
    M = 1; % assumed, should cancel out
    g = 386;
    wn = 2*pi/Tn(i);
    C = 2*E*wn*M;
    K = wn^2*M;
    uy = Cy(i)*M*g/K;
    [~,~,~,Sd(i),~,Sa(i)] = NewmarkAverageAcceleration(Tn(i), E, A_Parking, dt, u0, v0, uy);
    disp(['Tn = ' num2str(Tn(i))]);
    K*Sd(i)/(M*g)
end

%%

load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
u0 = 0;
v0 = 0;
Cy = [0.0859 0.0429];
Tn = [0.4 1.0];
E = 0.05;
for i = 1:2
    [u{i}, v{i}, a{i}, Sd(i), Sv(i), Sa(i), ~, ~, ~, mu(i)] = NewmarkAverageAccelerationCy(Tn(i), E, A_Parking, dt, u0, v0, Cy(i));
end


%% Part c
load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
Tn = 1;
E = 0.05;
u0 = 0;
v0 = 0;
fn = @(u) FindMu(Tn, E, A_Parking, dt, u0, v0, u) - 4;
M = 1; % assumed, should cancel out
g = 386;
wn = 2*pi/Tn;
C = 2*E*wn*M;
K = wn^2*M;

Cy = K*fzero(fn,0.5)/(M*g);


Ce = 0.01:0.01:1;
uy = Ce*M*g/K;
n = numel(uy);
muvec = zeros(n,1);
for i = 1:n
    muvec(i) = FindMu(Tn, E, A_Parking, dt, u0, v0, uy(i));
end

figure;
plot(muvec,Ce,'b-', [4 4], [min(Ce) max(Ce)], 'r--', 4, Cy, 'mo'); grid on;
xlabel('\mu'); ylabel('C_y');
xlim([0 6]); ylim([min(Ce) max(Ce)]);
title('8-Story T_n = 1.0s - Loma Prieta - Parking Record');
legend('Calculated \mu','\mu = 4', 'Calculated root');

%% Part e

load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
Tn = [0.4 1.0];
E = 0.05;
u0 = 0;
v0 = 0;
Ce = 0.01:0.01:1;
n = numel(Ce);
muvec = zeros(n,2);
for j = 1:2
    M = 1; % assumed, should cancel out
    g = 386;
    wn = 2*pi/Tn(j);
    C = 2*E*wn*M;
    K = wn^2*M;
    uy = Ce*M*g/K;
    for i = 1:n
        muvec(i,j) = FindMu(Tn(j), E, A_Parking, dt, u0, v0, uy(i));
    end
end

figure;
plot(muvec(:,1),Ce,'b-',muvec(:,2),Ce,'m--');
xlabel('\mu'); ylabel('C_y');
xlim([0 6]); ylim([min(Ce) max(Ce)]);
legend('T_n=0.4s','T_n=1.0s','Location','best');
title('Lateral Strength vs. Ductility');
grid on;

%% Part f

% Site Class C
a = 48;
b = 1.8;
c = 50;
Ts = 0.85;

% Selected R, T values
R = [2 3 4];
m = 100;
T = [0.4 1 2]%linspace(0.01,2,m);
m = 3;

% GM Information
load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
E = 0.05;
u0 = 0;
v0 = 0;

% Initialization
Sd = zeros(m,3);
Sd_inel = zeros(m,3);
CR = zeros(m,3);

for i = 1:3
    for j = 1:m
        CR(j,i) = 1 + ((1/(a*(T(j)/Ts)^b)) - 1/c)*(R(i)-1);
        [~,~,~,Sd(j,i),~,~] = NewmarkAverageAcceleration(T(j), E, A_Parking, dt, u0, v0, 1e10); % NOTE: uy sufficiently high so that remains elastic
        Sd_inel(j,i) = CR(j,i)*Sd(j,i);
    end
end

% Plotting
figure;
subplot(1,2,1); plot(T, Sd(:,1), 'r-', T, Sd(:,2), 'g-', T, Sd(:,3), 'b-'); grid on;
xlabel('Period, T [s]'); ylabel('Elastic Spectral Displacement [in]'); title('Linear Elastic Spectral Displacements');
legend('R=2','R=3','R=4','Location','best');
xlim([0 max(T)]); ylim([0 5.5]);
subplot(1,2,2); plot(T, Sd_inel(:,1), 'r-', T, Sd_inel(:,2), 'g-', T, Sd_inel(:,3), 'b-'); grid on;
xlabel('Period, T [s]'); ylabel('C_R*(Elastic Displacement) [in]'); title('Inelastic Spectral Displacements, Ruiz-Garcia & Miranda 2003');
legend('R=2','R=3','R=4','Location','best');
xlim([0 max(T)]); ylim([0 5.5]);

%% 

load Parking;
A_Parking = Parking(:,2)*0.00102;
dt = Parking(2,1);
u0 = 0;
v0 = 0;
Cy = [0.2 0.069];
Tn = [0.4 1.0];
E = 0.05;
for i = 1:2
    [u{i}, v{i}, a{i}, Sd(i), Sv(i), Sa(i), ~, ~, ~, mu(i)] = NewmarkAverageAccelerationCy(Tn(i), E, A_Parking, dt, u0, v0, Cy(i));
end
