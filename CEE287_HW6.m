%% CEE 287: Homework 6
% Nicole Paul

%% Part a

clear;
clc;

% Loading Data
load HW6Data;
A_rec = [A_3 A_7 A_R];

% Structure info
nfloors = 9;
E = 0.035;
floorofint = [3 7 10] - 1; % adjusting for labeling in HW
Ta = 2;
h_story = [18*12 11*12 12*12+6 12*12+6 12*12+6 12*12+6 12*12+6 12*12+6 13*12+7]';
h = (cumsum(h_story));
H = sum(h_story);

% Modal information - first 3 modes

% Period ratios
nmodes = 3;
T = zeros(1, nmodes);
for i = 1:nmodes
    T(i) = Ta/(2*i-1);
end

% Mode shapes
phi = zeros(nfloors, nmodes);
for i = 1:nmodes
    phi(:,i) = (-1)^(i-1)*sin(((2*i-1)*pi*h)./(2*H));
end

% Modal participation factors
Gamma = zeros(1, nmodes);
for i = 1:nmodes
    Gamma(i) = (4*(-1)^(i-1))/(2*(i*pi)-pi);
end

% Running (elastic) time history analysis with approximate periods
conv_to_g = 0.00102;
u0 = 0; v0 = 0; isPlot = 0; dt = t(2); n = numel(t);
ag_rel = zeros(n,nmodes); ag_abs = zeros(n, nmodes);
Sam = zeros(1, nmodes);
for i = 1:nmodes
    [~,~,a,~,~,Sa,~,~] = RecurrenceSDOF(T(i),E,A_g*conv_to_g,dt,u0,v0,isPlot);
    ag_rel(:,i) = a-A_g*conv_to_g;
    ag_abs(:,i) = a;
    Sam(i) = Sa;
end

% Using modal information to get history for each floor of interest
n_floorofint = numel(floorofint);
Am_floor_rel = zeros(n,nmodes,n_floorofint);
A_floor_rel = zeros(n, n_floorofint);
Am_floor_abs = zeros(n,nmodes,n_floorofint);
A_floor_abs = zeros(n, n_floorofint);
for j = 1:n_floorofint
    f = floorofint(j);
    for i = 1:nmodes
        Am_floor_rel(:,i,j) = Gamma(i)*phi(f,i)*ag_rel(:,i);
        Am_floor_abs(:,i,j) = Gamma(i)*phi(f,i)*ag_abs(:,i);
    end
    A_floor_rel(:,j) = A_g*conv_to_g + sum(Am_floor_rel(:,:,j),2);
    A_floor_abs(:,j) = (1-Gamma*phi(f,:)')*A_g*conv_to_g + sum(Am_floor_abs(:,:,j),2);
end

% Generating figures
for j = 1:n_floorofint
    figure;
    plot( t, conv_to_g*A_rec(:,j),'k-', t, A_floor_rel(:,j), 'b-', t, A_floor_abs(:,j), 'g--' ,'LineWidth', 1);
    xlabel('Time, t [s]'); ylabel('Acceleration [g]'); title(['Floor ' num2str(floorofint(j)+1) ' Acceleration Time History'] );
    legend('Recorded','Approximate Relative','Approximate Absolute');
    grid on;
    xlim([0 t(end)]); ylim([-0.2 0.2]);
end

%% Part b

% Using modal information to get history for each floor of interest
floorofint = 1:9; n_floorofint = numel(floorofint);
Am_floor_rel = zeros(n,nmodes,n_floorofint); A_floor_rel = zeros(n, n_floorofint);
Am_floor_abs = zeros(n,nmodes,n_floorofint); A_floor_abs = zeros(n, n_floorofint);
Am_floor_rec = zeros(n_floorofint, nmodes); A_floor_rec = zeros(n_floorofint,1);
PFA_rel = zeros(nfloors, 1); PFA_abs = zeros(nfloors, 1); PFA_rec = zeros(nfloors, 1);
for j = 1:n_floorofint
    for i = 1:nmodes
        Am_floor_rel(:,i,j) = Gamma(i)*phi(j,i)*ag_rel(:,i);
        Am_floor_abs(:,i,j) = Gamma(i)*phi(j,i)*ag_abs(:,i);
        Am_floor_rec(j,i) = (Gamma(i)*phi(j,i)*Sam(i))^2;
    end
    A_floor_rel(:,j) = A_g*conv_to_g + sum(Am_floor_rel(:,:,j),2);
    A_floor_abs(:,j) = (1-Gamma*phi(j,:)')*A_g*conv_to_g + sum(Am_floor_abs(:,:,j),2);
    A_floor_rec(j) = sqrt(((1-Gamma*phi(j,:)')*max(abs(A_g))*conv_to_g)^2 + sum(Am_floor_rec(j,:)));
    
    PFA_rel(j) = max(abs(A_floor_rel(:,j)));
    PFA_abs(j) = max(abs(A_floor_abs(:,j)));
    PFA_rec(j) = A_floor_rec(j);
end
PGA = max(abs(A_g))*conv_to_g;

% ASCE7-10 calc
ap = 1;
Rp = 1;
I = 1;
PFA_code = PGA*(ap/(Rp/I))*(1+2*h/H);

% Plotting results
floors = 1:10;
figure;
plot([PGA; PFA_rel], floors, 'bo-', [PGA; PFA_abs], floors, 'go--', [PGA; PFA_rec], floors, 'ko-', [PGA; PFA_code], floors, 'ro-','LineWidth', 1); grid on;
xlabel('Peak Floor Acceleration, PFA [g]'); ylabel('Floor');
title('Peak Floor Accelerations');
legend('RHA Relative', 'RHA Absolute', 'RSA Absolute', 'ASCE7-10', 'Location', 'best');
xlim([0 0.3]); ylim([1 10]);

%% Part c

% Floors of interest
floorofint = [3 7 10] - 1; % adjusting for labeling in HW
n_floorofint = numel(floorofint);

% Range of periods for response spectrum
m = 1000;
Trange = linspace(0,4,m);

% Running (elastic) time history analysis with approximate periods
FloorSa_rel = zeros(m,n_floorofint);
FloorSa_abs = zeros(m,n_floorofint);
FloorSa_rec = zeros(m,n_floorofint);
u0 = 0; v0 = 0; isPlot = 0; dt = t(2); n = numel(t);
ag_rel = zeros(n,nmodes); ag_abs = zeros(n, nmodes);
E = 0.05;
for i = 1:n_floorofint
    for j = 1:m
        f = floorofint(i);
        % Using approx absolute
        [~,~,~,~,~,Sa,~,~] = RecurrenceSDOF(Trange(j),E,A_floor_rel(:,f),dt,u0,v0,isPlot);
        FloorSa_rel(j,i) = Sa;
        % Using approx relative
        [~,~,~,~,~,Sa,~,~] = RecurrenceSDOF(Trange(j),E,A_floor_abs(:,f),dt,u0,v0,isPlot);
        FloorSa_abs(j,i) = Sa;
        % Using recording
        [~,~,~,~,~,Sa,~,~] = RecurrenceSDOF(Trange(j),E,A_rec(:,i)*conv_to_g,dt,u0,v0,isPlot);
        FloorSa_rec(j,i) = Sa;
    end
end

figure;
% Generating figures
for j = 1:n_floorofint
    subplot(1,3,j);
    plot(Trange, FloorSa_rec(:,j),'k-', Trange, FloorSa_rel(:,j), 'b-', Trange, FloorSa_abs(:,j), 'g--' ,'LineWidth', 1);
    xlabel('Period, T [s]'); ylabel('Spectral Acceleration, S_a [g]'); title(['Floor ' num2str(floorofint(j)+1) ' Spectra'] );
    legend('Recorded','Relative','Absolute','Location','NorthWest');
    grid on;
    xlim([0 4]); ylim([0 1.2]);
end

% Generating figures
for j = 1:n_floorofint
    subplot(1,3,j);
    f = floorofint(j);
    plot(Trange, FloorSa_rec(:,j)./PFA_rec(f),'k-', Trange, FloorSa_rel(:,j)./PFA_rel(f), 'b-', Trange, FloorSa_abs(:,j)./PFA_abs(f), 'g--' ,[0 4], [2.5 2.5],'r--');
    xlabel('Period, T [s]'); ylabel('Normalized Spectral Acceleration, S_a/PFA'); title(['Normalized Floor ' num2str(floorofint(j)+1) ' Spectra'] );
    legend('Recorded','Relative','Absolute','ASCE','Location','NorthWest');
    grid on;
    xlim([0 4]); ylim([0 8]);
end
% % Plotting
% figure;
% plot(Trange,FloorSa); grid on;
% xlabel('Period, T [s]'); ylabel('Spectral Acceleration, S_a [g]');
% title('Elastic Response Spectra for Floor Acceleration History for\xi = 5%');
% legend('Floor 3', 'Floor 7', 'Roof', 'Location', 'best');