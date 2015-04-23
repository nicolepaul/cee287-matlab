%% CEE 287: Homework 3
% Nicole Paul

%% Problem 1

clear;
clc;

% From USGS, for site
Ss = 1.767;
S1 = 0.798;
TL = 12;

% Table 11.4-1 and Table 11.4-2
Fa = 1.0;
Fv = 1.5;

% Equations 11.4-1 and 11.4-2
Sms = Fa*Ss;
Sm1 = Fv*S1;

% Equations 11.4-3 and 11.4-4
Sds = 2*Sms/3;
Sd1 = 2*Sm1/3;

% Section 11.4.5 Design Response Spectrum
T0 = 0.2*Sd1/Sds;
Ts = Sd1/Sds;

% Calculating Design Spectra
m = 1000;
T = linspace(0,14,m)';
Sa = zeros(m,1);
for i = 1:m
    if T(i) < T0
        Sa(i) = Sds*(0.4+0.6*T(i)/T0);
    elseif T(i) >= T0 && T(i) < Ts
        Sa(i) = Sds;
    elseif T(i) >= Ts && T(i) < TL
        Sa(i) = Sd1/T(i);
    else
        Sa(i) = Sd1*TL/T(i)^2;
    end
end

% Plotting
figure;
plot(T,Sa,'b-',T,1.5*Sa,'r--'); grid on;
xlabel('Period, T [s]');
ylabel('Spectral Acceleration, S_a [g]');
xlim([0 2]);
legend('DBE','MCE','Location','best');
title('Response Spectrum');
        
%% Problem 3

% Getting ground motion information
load('Parking.mat');
load('SLAC_1.mat');
load('SLAC_2.mat');
load('VA_1.mat');
load('VA_2.mat');

convtog = 0.00101971621;

GM = convtog*[Parking(:,2) VA_1(:,2) VA_2(:,2) SLAC_1(:,2) SLAC_2(:,2)];
t = Parking(:,1); % same for all GM
nGM = 5;

% Period Range
m = 1000;
T_range = linspace(0,4,m);

% Initial conditions
u0 = 0;
v0 = 0;

% Other Inputs
dt = t(2);
E = 0.05;

% Initialization
Sd = zeros(m,nGM);
Sa = zeros(m,nGM);

% Getting response spectrum for each GM
for i = 1:nGM
    for j = 1:m
        [~,~,~,Sd(j,i),~,Sa(j,i),~,~] = RecurrenceSDOF(T_range(j),E,GM(:,i),dt,u0,v0);
    end
end

% Range of periods of interest
Tn = 0.6;
T_lower = 0.2*Tn;
T_upper = 1.5*Tn;

% % % DESIGN SPECTRA
% From USGS, for site
Ss = 1.767;
S1 = 0.798;
TL = 12;

% Table 11.4-1 and Table 11.4-2
Fa = 1.0;
Fv = 1.5;

% Equations 11.4-1 and 11.4-2
Sms = Fa*Ss;
Sm1 = Fv*S1;

% Equations 11.4-3 and 11.4-4
Sds = 2*Sms/3;
Sd1 = 2*Sm1/3;

% Section 11.4.5 Design Response Spectrum
T0 = 0.2*Sd1/Sds;
Ts = Sd1/Sds;
inds_T0 = find(T_range<T_lower);
inds_Ts = find(T_range>T_upper);
ind_T0 = inds_T0(end);
ind_Ts = inds_Ts(1);

% Calculating Design Spectra
Sa_DBE = zeros(m,1);
for i = 1:m
    if T_range(i) < T0
        Sa_DBE(i) = Sds*(0.4+0.6*T_range(i)/T0);
    elseif T_range(i) >= T0 && T_range(i) < Ts
        Sa_DBE(i) = Sds;
    elseif T_range(i) >= Ts && T_range(i) < TL
        Sa_DBE(i) = Sd1/T_range(i);
    else
        Sa_DBE(i) = Sd1*TL/T_range(i)^2;
    end
end

% Finding appropriate scale factors
range = ind_T0:ind_Ts;
SF = Sa_DBE(range)./mean(Sa(range,:),2);
max_SF = max(SF);

% Plotting
figure;
plot(T_range, Sa*max_SF); grid on; hold on;
plot(T_range, mean(Sa,2)*max_SF, 'k--', 'LineWidth', 2);
plot(T_range, Sa_DBE, 'r--', 'LineWidth', 2);
xlabel('Period, T [s]');
ylabel('Spectral Acceleration, S_a [g]');
title(['Response Spectrum, Scale Factor = ' num2str(max_SF)]);
legend('Parking','VA_1','VA_2','SLAC_1','SLAC_2','Mean','DBE');


%% Problem 4
%% Problem 4

% Getting ground motion information
load('Parking.mat');
load('SLAC_1.mat');
load('SLAC_2.mat');
load('VA_1.mat');
load('VA_2.mat');

convtog = 0.00101971621;

GM = convtog*[Parking(:,2) VA_1(:,2) VA_2(:,2) SLAC_1(:,2) SLAC_2(:,2)];
t = Parking(:,1); % same for all GM
nGM = 5;
%%
% Period Range
m = 1000;
T_range = linspace(0,6,m);

% Initial conditions
u0 = 0;
v0 = 0;

% Other Inputs
dt = Parking(2,1);
E = 0.05;

% Initialization
Sd = zeros(m,nGM);
Sa = zeros(m,nGM);

% Getting response spectrum for each GM
for i = 1:nGM
    for j = 1:m
        [~,~,~,Sd(j,i),~,Sa(j,i),~,~] = RecurrenceSDOF(T_range(j),E,GM(:,i),dt,u0,v0);
    end
end

% Range of periods of interest
Tn = 0.6;
T_lower = 0.2*Tn;
T_upper = 1.5*Tn;

% % % DESIGN SPECTRA
% From USGS, for site
Ss = 1.767;
S1 = 0.798;
TL = 12;

% Table 11.4-1 and Table 11.4-2
Fa = 1.0;
Fv = 1.5;

% Equations 11.4-1 and 11.4-2
Sms = Fa*Ss;
Sm1 = Fv*S1;

% Equations 11.4-3 and 11.4-4
Sds = 2*Sms/3;
Sd1 = 2*Sm1/3;

% Section 11.4.5 Design Response Spectrum
T0 = 0.2*Sd1/Sds;
Ts = Sd1/Sds;
inds_T0 = find(T_range<T_lower);
inds_Ts = find(T_range>T_upper);
ind_T0 = inds_T0(end)+1;
ind_Ts = inds_Ts(1)-1;

% Calculating Design Spectra
Sa_DBE = zeros(m,1);
for i = 1:m
    if T_range(i) < T0
        Sa_DBE(i) = Sds*(0.4+0.6*T_range(i)/T0);
    elseif T_range(i) >= T0 && T_range(i) < Ts
        Sa_DBE(i) = Sds;
    elseif T_range(i) >= Ts && T_range(i) < TL
        Sa_DBE(i) = Sd1/T_range(i);
    else
        Sa_DBE(i) = Sd1*TL/T_range(i)^2;
    end
end

% Finding appropriate scale factors
range = ind_T0:ind_Ts;
SF = 1.5*Sa_DBE(range)./mean(Sa(range,:),2);
max_SF = max(SF);
%%
% Plotting
figure;
plot(T_range, Sa*max_SF); grid on; hold on;
plot(T_range, mean(Sa,2)*max_SF, 'k--', 'LineWidth', 2);
plot(T_range, 1.5*Sa_DBE, 'r--', 'LineWidth', 2);
xlabel('Period, T [s]');
ylabel('Spectral Acceleration, S_a [g]');
title(['Response Spectrum, Scale Factor = ' num2str(max_SF)]);
legend('Parking','VA_1','VA_2','SLAC_1','SLAC_2','Mean','MCE');


%% Problem 5


Tn = 0.6;
E = 0.02;

SF_DBE = 2.9301;
SF_MCE = 4.3852;

Sd_DBE = zeros(nGM,1);
Sa_DBE = zeros(nGM,1);
Sd_MCE = zeros(nGM,1);
Sa_MCE = zeros(nGM,1);
for i = 1:nGM
    [~,~,~,Sd_DBE(i),~,Sa_DBE(i),~,~] = RecurrenceSDOF(Tn,E,GM(:,i)*SF_DBE,dt,u0,v0);
    [~,~,~,Sd_MCE(i),~,Sa_MCE(i),~,~] = RecurrenceSDOF(Tn,E,GM(:,i)*SF_MCE,dt,u0,v0);
end

R = 8;
Cd = 5;

del_DBE = (Cd/R)*Sd_DBE
del_MCE = (Cd/R)*Sd_MCE

%% Calculating R values

R = 8;
Tn = 0.6;
Fy = 1.178/R;
% Cy = 0.1841;
% Fy = Cy;

R_DBE = Sa_DBE/Fy
R_MCE = Sa_MCE/Fy


%% Problem 6

% Site Class D
a = 57;
b = 1.85;
c = 60;
Ts = 1.05;

% Selected R, T values
R = 8;
T = 0.6;

% GM Information
dt = Parking(2,1);
E = 0.05;
u0 = 0;
v0 = 0;

% Initialization
CR = zeros(nGM,1);
Sd_DBE = zeros(nGM,1);
Sd_inel_DBE = zeros(nGM,1);
Sd_MCE = zeros(nGM,1);
Sd_inel_MCE = zeros(nGM,1);
for i = 1:nGM
    CR_DBE(i) = 1 + ((1/(a*(T/Ts)^b)) - 1/c)*(R_DBE(i)-1);
    CR_MCE(i) = 1 + ((1/(a*(T/Ts)^b)) - 1/c)*(R_MCE(i)-1);
    [~,~,~,Sd_DBE(i),~,~] = NewmarkAverageAcceleration(T, E, SF_DBE*GM(:,i), dt, u0, v0, 1e10); % NOTE: uy sufficiently high so that remains elastic
    Sd_inel_DBE(i) = CR_DBE(i)*Sd_DBE(i);
    [~,~,~,Sd_MCE(i),~,~] = NewmarkAverageAcceleration(T, E, SF_MCE*GM(:,i), dt, u0, v0, 1e10); % NOTE: uy sufficiently high so that remains elastic
    Sd_inel_MCE(i) = CR_MCE(i)*Sd_MCE(i);
end

%% Problem 7

Tn = 0.6;
m = 1;
k = m*(2*pi)^2/Tn^2;
R = 8;
Fy = 1.178/R;
uy = Fy./k;
E = 0.02;

for i = 1:nGM
    [~,~,~,Sdi_DBE(i),~,~] = NewmarkAverageAcceleration(Tn, E, SF_DBE*GM(:,i), dt, u0, v0, uy); % NOTE: uy sufficiently high so that remains elastic
    [~,~,~,Sdi_MCE(i),~,~] = NewmarkAverageAcceleration(Tn, E, SF_MCE*GM(:,i), dt, u0, v0, uy); % NOTE: uy sufficiently high so that remains elastic
end

%% Problem 7

Tn = 0.6;
m = 1;
k = m*(2*pi)^2/Tn^2;
R = 8;
Fy = 1.178/R;
uy = Fy./k;
E = 0.02;
Cy = 0.1841;

for i = 1:nGM
    [~,~,~,Sdi_DBE(i),~,~] = NewmarkAverageAccelerationCy(Tn, E, SF_DBE*GM(:,i), dt, u0, v0, Cy); % NOTE: uy sufficiently high so that remains elastic
    [~,~,~,Sdi_MCE(i),~,~] = NewmarkAverageAccelerationCy(Tn, E, SF_MCE*GM(:,i), dt, u0, v0, Cy); % NOTE: uy sufficiently high so that remains elastic
end
