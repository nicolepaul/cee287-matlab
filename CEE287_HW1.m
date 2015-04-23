%% CEE 287: Homework 1
% Max Ferguson
% Nicole Paul

%% PART A

%% Problem 1

load ElCentroGMTH
Tn = 1;
E = 0.02;
u0 = 0;
v0 = 0;
[u,v,at,Sd,Sv,Sa,PSv,PSa] = RecurrenceSDOF(Tn,E,A_ElCentro,dt,u0,v0);

%% Problem 2

E = [0.02 0.05];
Tn = [0.3 1.0 2.0];
dt = 0.005;
conv_g = 0.00102;

load VA_2
A_V2 = VA_2(:,2)*conv_g;

u0 = 0;
v0 = 0;

t = 0:dt:dt*numel(A_V2)-dt;
for i = 1:numel(Tn)
    for j = 1:numel(E)
        [u,v,a,Sd,Sv,Sa,PSv,PSa] = RecurrenceSDOF(Tn(i),E(j),A_V2,dt,u0,v0);
        t_Sd = find(abs(u)==Sd);
        t_Sv = find(abs(v)==Sv);
        t_Sa = find(abs(a)==Sa);
        ta_peaks{i,j} = [Sd Sv Sa];
        ta_times{i,j} = [t(t_Sd) t(t_Sv) t(t_Sa)];
    end
end

%% Problem 5

load BJFData
load LPRecurrenceData

b1 = B1SS;
b2 = B2;
b3 = B3;
b5 = B5;
bv = BV;
Vs = 1116.2*0.3048; % m/s
Rjb = 28; % km
M = 6.9;

nT = numel(Period);
lnY = zeros(nT,1);
Y = zeros(nT,1);
Y_15 = zeros(nT,1);
Y_85 = zeros(nT,1);
for i = 1:nT
    r = sqrt(Rjb^2 + h(i)^2);
    lnY(i) = b1(i) + b2(i)*(M-6) + b3(i)*(M-6)^2 + b5(i)*log(r) + bv(i)*log(Vs/VA(i));
    Y(i) = exp(lnY(i));
    Y_15(i) = exp(lnY(i) - (slnY(i)));
    Y_85(i) = exp(lnY(i) + (slnY(i)));
end

figure;
plot(Period,Y,'ko-',Period,Y_15,'kv-',Period,Y_85,'k^-',Period,handles.PSa); grid on; hold on;
plot(Period,mean(handles.PSa,2),'k--','LineWidth',2); ylim([0 1.2]);
xlabel('Period, T [s]');
ylabel('Spectral Acceleration, S_a [g]');
title('Response Spectrum as per BJF 1997');
legend('GMPE 50%','GMPE 15%','GMPE 85%','Parking','SLAC_1','SLAC_2','VA_1','VA_2','Average of GM','Location','best');


%% Problem 6

load LPRecurrenceData

load BJFData

b1 = B1SS;
b2 = B2;
b3 = B3;
b5 = B5;
bv = BV;
Vs = 1116.2*0.3048; % m/s
Rjb = 28; % km
M = 6.9;

nT = numel(Period);
lnY = zeros(nT,1);
Y = zeros(nT,1);
nGM = size(handles.PSa,2);
eps_Y = zeros(nT,nGM);
for i = 1:nT
    r = sqrt(Rjb^2 + h(i)^2);
    lnY(i) = b1(i) + b2(i)*(M-6) + b3(i)*(M-6)^2 + b5(i)*log(r) + bv(i)*log(Vs/VA(i));
    Y(i) = exp(lnY(i));
    for j = 1:nGM % for each ground motion
        eps_Y(i,j) = (log(handles.PSa(i,j)) - lnY(i))/slnY(i);
    end
end

figure;
plot(Period,eps_Y(:,1),'r*-',Period,eps_Y(:,2),'m*-',Period,eps_Y(:,3),'g*-',Period,eps_Y(:,4),'b*-',Period,eps_Y(:,5),'c*-'); grid on; hold on;
plot(Period,mean(eps_Y,2),'k--','LineWidth',2);
legend('Parking','SLAC_1','SLAC_2','VA_1','VA_2','Mean','Location','best');
xlabel('Period, T [s]');
ylabel('Epsilon');
title('Epsilon Spectrum for Loma Prieta, BJF 1997');
ylim([0 3]); xlim([0 2.2]);

figure;
plot(Period, slnY./abs(lnY),'ro-',Period, slnY,'mo-',handles.Tn, sig_lnY,'b*-'); grid on;
legend('\sigma_{ln}/\mu_{ln} of GMPE','\sigma_{ln} of GMPE', 'Variability of 5 GM');
title('Comparison of Dispersion of GMPE and 5 GM'); ylim([0 1]);
xlabel('Period, T [s]'); ylabel('Uncertainty');
