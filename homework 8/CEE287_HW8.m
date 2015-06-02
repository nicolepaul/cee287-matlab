%% CEE 287: Homework 8
% Max Ferguson, Nicole Paul

%% Part A

%syms alpha real positive
close all;
alpha =  0.056409983382519;

gravity = 386.1; % in/s^2
designSd1 = 1.0;
targetDrift = 0.5/100; % 0.5%

% Define the height of each floor in the building
% The vector is defined from the bottom to top, starting at 0'
nfloors = 9;
hroof = 118;
hfloor = 118/9;
heights = transpose(0:hfloor:hroof);

% Set the building parameters as we they were defined
% in homework 5, and redefined in homework 8.
W = 1800; % kips/in
mass = W/386.1; % At each floor
stiffness = 1700; % At each floor
nfloors = 9;
nmodes = 1;

% Get the stiffness and mass matrix of the equivalent structure
% Compute the mode shape for the first mode of vibration
[M, K] = computeMatrices(nfloors,mass,stiffness);
[~,T,sphi,GammaStruct] = eigenvalueAnalysis(nfloors,nmodes,mass,stiffness);

% Use the modal analysis to compute the equivalent lateral stiffness
% and mass of the structure. These parameters will be normalization
% dependant but the ratio of them will be unique.
beta = 1/nfloors; % Assume same mass at each floor

fprintf('Part A:\n')
ks = sphi(:,1)'*K*sphi(:,1);
ms = sphi(:,1)'*M*sphi(:,1);
ki = alpha*ks;
mi = beta*(nfloors*mass);

fprintf('Equivalent Stiffness, ks = %.4f [kips/in per unit mass]\n',ks);
fprintf('Equivalent Mass, ms = %.4f [kips.s^2/in per unit mass]\n',ms);
% fprintf('Isolator Stiffness, ki = %.4f [kips/in per unit mass]\n',ki);
% fprintf('Isolator Mass, mi = %.4f [kips.s^2/in per unit mass]\n',mi);

% Analyze the equivalent 2DOF structure with the following properties
% ks,ms,ki,mi. First a modal analysis is conducted to find the period
% of each mode. Then another modal analysis is conducted to find the
% displacement and drift at each floor.
w1squared = (ms*ks*(1+alpha)+beta*ms*ks - sqrt((ms*ks*(1+alpha)+beta*ms*ks)^2 - 4*alpha*beta*ms^2*ks^2))/(2*beta*ms^2);
w2squared = (ms*ks*(1+alpha)+beta*ms*ks + sqrt((ms*ks*(1+alpha)+beta*ms*ks)^2 - 4*alpha*beta*ms^2*ks^2))/(2*beta*ms^2);

T1 = 2*pi/sqrt(w1squared);
T2 = 2*pi/sqrt(w2squared);

gamma1 = (-ms*w1squared+ks)/ks;
gamma2 = (-ms*w2squared+ks)/ks;

phi1 = [1, gamma1]';
phi2 = [1, gamma2]';

Gamma1 = (1+gamma1*beta)/(1+gamma1^2*beta);
Gamma2 = (1+gamma2*beta)/(1+gamma2^2*beta);

E1 = 0.1;
E2 = 0.035;

An1 = linChangDamping(designSd1/T1, T1, E1);
An2 = linChangDamping(designSd1, T2, E2);

Uj1 = Gamma1*phi1*An1*gravity/w1squared;
Uj2 = Gamma2*phi2*An2*gravity/w2squared;

Beta1 = GammaStruct*sphi(end);
temp_sphi = [0; sphi];
Beta2 = max(hroof*(temp_sphi(2:end)-temp_sphi(1:end-1))./(heights(2:end)'*sphi(end))');

u = sqrt(Uj1.^2+Uj2.^2);
Sd1 = abs(u(2) - u(1));
IDR_maxSDOF = Sd1/((2/3)*hroof*12);
IDR_maxMDOF = (2/3)*Beta1*Beta2*IDR_maxSDOF;

fprintf('The drift %.4f [%%]\n',100*IDR_maxMDOF)
%alpha_design = eval(solve(IDR_maxMDOF-targetDrift==0))

%% 

% targetdrift = 0.005;
% f = @(alpha) targetdrift - get_drift(alpha);
% alpha_design = fsolve2(f, 0.05);
% 
% m = 1000;
% alpha_range = linspace(0,2,m);
% drift_range = NaN(size(alpha_range));
% for i = 1:m
%     drift_range(i) = get_drift(alpha_range(i));
% end
% 
% figure;
% plot(alpha_range, drift_range, 'b-', alpha_design, targetdrift, 'ro'); grid on;
% xlabel('Alpha, \alpha'); ylabel('Max Interstory Drift Ratio');
% title('Influence of Alpha on IDR_{max}');


%% Part B

% Scale up stiffness from 2DOF system
MDOF_mass = nfloors*mass;
MDOF_stiffness = 16.6; %%%%%%%%%%%ks * MDOF_mass/ms;
stiffness_fixed = stiffness*ones(1,nfloors);
stiffness_isolated = [MDOF_stiffness,stiffness_fixed];

mass_fixed = mass*ones(1,nfloors);
mass_isolated = mass*ones(1,nfloors+1);


% 1) Get the stiffness and mass matrix of the equivalent structure
[M_fixed, K_fixed] = computeMatrices(nfloors, mass_fixed, stiffness_fixed);
[M_isolated, K_isolated] = computeMatrices(nfloors+1, mass_isolated, stiffness_isolated);

% Determining fundamental periods
[~,T_fixed,sphi_fixed,Gamma_fixed] = eigenvalueAnalysis(nfloors,nfloors,mass_fixed,stiffness_fixed);
[~,T_isolated,sphi_isolated,Gamma_isolated] = eigenvalueAnalysis(nfloors+1,nfloors+1,mass_isolated,stiffness_isolated);

% Determining C_sm (seismic coefficient for each mode)
Sds = 1.5;
Sd1 = 1.0;
A_fixed = min([Sds*ones(nfloors,1) Sd1./T_fixed],[],2);
A_isolated = min([Sds*ones(nfloors+1,1) Sd1./T_isolated],[],2);

% Modify for damping ratios
Ef = [E2,E2*ones(1,8)];
Eb = [E1,E2*ones(1,9)];
Csm_fixed = linChangDamping(A_fixed,T_fixed',Ef);
Csm_isolated = linChangDamping(A_isolated,T_isolated',Eb);

% Performing response spectrum analysis
[~,V_fixed,U_fixed,drift_fixed] = modalAnalysis(nfloors,mass_fixed,stiffness_fixed,Csm_fixed,hfloor);
[~,V_isolated,U_isolated,drift_isolated] = modalAnalysis(nfloors+1,mass_isolated,stiffness_isolated,Csm_isolated,hfloor);

% Calculate total acceleration according to Miranda,
% Estimation of acceleration demands in structures, Page 22.
% As RecuranceSDOF provides absolute acceleration we will use method 2
term1 = zeros(nfloors,length(A_fixed));
term2 = zeros(nfloors,length(A_fixed));
PGA = 0.4*Sds;
for i=1:length(A_fixed)
    GammaPhi = Gamma_fixed(i)*sphi_fixed(:,i);
    term1(:,i) = GammaPhi;
    term2(:,i) = (GammaPhi*A_fixed(i)).^2;
end
A_fixed_absolute = sqrt( ((1-sum(term1,2))*PGA).^2 + sum(term2,2) );

% Calculate total acceleration according to Miranda,
% Estimation of acceleration demands in structures, Page 22.
% As RecuranceSDOF provides absolute acceleration we will use method 2
term1 = zeros(nfloors,length(A_isolated));
term2 = zeros(nfloors,length(A_isolated));
for i=1:length(A_isolated)
    GammaPhi = Gamma_isolated(i)*sphi_isolated(:,i);
    term1(:,i) = GammaPhi;
    term2(:,i) = (GammaPhi*Ai(i)).^2;
end
A_isolated_absolute = sqrt( ((1-sum(term1,2))*PGA).^2 + sum(term2,2) );





V_fixed = srss(V_fixed);
U_fixed = vertcat(srss(U_fixed),0);
drift_fixed = srss(drift_fixed);
V_isolated = srss(V_isolated);
U_isolated = srss(U_isolated);
drift_isolated = srss(drift_isolated);

%  2) Comparison plot of Story displacements for isolated vs fixed base structure
figure; hold on;
plot(U_fixed,flip(heights),'-o');
plot(U_isolated,flip(heights),'-o');
legend({'Fixed','Isolated'},'Location','best');
title('Displacement')
xlabel('Displacement [in]');
ylabel('Height [ft]');

% 2a) Comparison plot of superstructure displacement
figure; hold on;
plot(U_fixed,flip(heights),'-o');
plot(U_isolated-U_isolated(end),flip(heights),'-o');
legend({'Fixed','Isolated'},'Location','best');
title('Superstructure Displacement')
xlabel('Displacement [in]');
ylabel('Height [ft]');

% 3) Plot the Story drift
% Lets assume the isolator is 12" high and recalculate the drift
drift_isolated(end) = (U_isolated(end-1)-U_isolated(end))/12;
figure; hold on;
plotSquare(100*drift_fixed,flip(heights(2:end)),'-o');
plotSquare(100*drift_isolated,flip(heights),'-o');
plot(0.5*ones(size(heights)),heights,':r')
legend({'Fixed','Isolated'},'Location','best');
title('Interstory Drift')
xlabel('Drift [%]');
ylabel('Height [ft]');
xlim([0,3.0]);


% 3) Plot the Story shear forces
% Lets assume the isolator is 12" high and recalculate the drift
figure; hold on;
stories = 1:length(heights)-1;
plot(V_fixed,flip(stories),'-o');
plot(V_isolated,flip([1 stories]),'-o');
legend({'Fixed','Isolated'},'Location','best');
title('Story Shear')
xlabel('Story shear [kips]');
ylabel('Story');

% Plot the Story shear forces
% Lets assume the isolator is 12" high and recalculate the drift
figure(5); hold on;
stories = 1:length(heights)-1;
plot(A_fixed, stories, '-o');
plot(A_isolated, [1 stories],'-o');
legend({'RSA Fixed','RSA Isolated','RSA Fixed','RSA Isolated'},'Location','best');
title('Story Shear')
xlabel('Floor acceleration [g]');
ylabel('Story');

% 3) Calculate the absolute floor acceleration using Miranda and Taghavi(2005) 
% This method is also described in the class notes, estimation of acceleration
% demands in structures, Page 22.

%term1 = zeros(nfloors,length(Csm_fixed));
%term2 = zeros(nfloors,length(Csm_fixed));
%for i=1:length(Ai)
%    GammaPhi = Gamma_fixed(i)*sphi_fixed(:,i);
%    term1(:,i) = GammaPhi;
%    term2(:,i) = (GammaPhi*Csm_fixed(i)).^2;
%end
%acceleration = sqrt( ((1-sum(term1,2))*max(A)).^2 + sum(term2,2) );
%plot([PGA acceleration'],floor,'ko-')





%% Part C
% Conduct  a  modal  response  history  analysis  on  both  the  isolated  and 
% fixed  base  structure  using  the  same ground motion from Assignment 6,
% but now scaled by a factor of 5.  Report the following:
scale_gm = 5;
conv_cm_to_g = 0.001019716;
hw6acc = xlsread('hw6_acc.xlsx');
hw6acc(:,2:5) = scale_gm*conv_cm_to_g*hw6acc(:,2:5);
    
A = hw6acc(:,2);
dt = hw6acc(2,1)-hw6acc(1,1);
    
% Modal response using stiffness matrices %%%%%
[M_isolated, ~] = computeMatrices(nfloors+1,mass,stiffness);
[uf,~,af,ff,vf,~,~,~] = modalTimeHistoryAbsolute(M,T_fixed,Gamma_fixed,sphi_fixed,Ef,A,dt);
[ub,~,ab,fb,vb,~,~,~] = modalTimeHistoryAbsolute(M_isolated,T_isolated,Gamma_isolated,sphi_isolated,Eb,A,dt);

% Calculate the drift
% This subtracts the i+1 floor displacement from the ith floor displacement
drift_fixed = uf-[zeros(length(uf),1) uf(:,1:end-1)];
drift_isolated = ub-[zeros(length(ub),1) ub(:,1:end-1)];

for i=1:length(uf)
    displacement = uf(i,:);
    drift_fixed(i,:) = displacement - [0 displacement(1:end-1)];
end

for i=1:length(ub)
    displacement = ub(i,:);
    drift_isolated(i,:) = displacement - [0 displacement(1:end-1)];
end

% Calculate the response envelopes
actual(1) = max(abs(hw6acc(:,2)));
F_fixed = max(abs(ff))';
V_fixed = max(abs(vf))';
A_fixed = max(abs(af))';
U_fixed = max(abs(uf))';
U_fixed = vertcat(0,U_fixed);
Drift_fixed = max(abs(drift_fixed))'/(12*hfloor);

F_isolated = max(abs(vf))';
V_isolated = max(abs(vb))';
A_isolated = max(abs(ab))';
U_isolated = max(abs(ub))';
U_isolator = ub(:,1);
U_normalized = max(ub-diag(U_isolator)*ones(size(ub)));
Drift_isolated = max(abs(drift_isolated))'/(12*hfloor);


% RHA: Comparison plot of Story displacements for isolated vs fixed base structure
figure(1); hold on;
plot(U_fixed,heights,'-o');
plot(U_isolated,heights,'-o');
title('Displacement')
xlabel('Displacement [in]');
ylabel('Height [ft]');
legend({'RSA Fixed','RSA Isolated','RSA Fixed','RSA Isolated'},'Location','best');


% RHA: Comparison plot of Story displacements for isolated vs fixed base structure
figure(2); hold on;
plot(U_fixed,heights,'-o');
plot(U_normalized,heights,'-o');
legend({'Fixed','Isolated'},'Location','best');
xlabel('Displacement [in]');
ylabel('Height [ft]');
legend({'RSA Fixed','RSA Isolated','RSA Fixed','RSA Isolated'},'Location','best');


% RHA: Comparison plot of superstructure displacement
figure(3); hold on;
plotSquare(100*flip(Drift_fixed), flip(heights(2:end)),'-o');
plotSquare(100*flip(Drift_isolated), flip(heights),'-o');
legend({'RSA Fixed','RSA Isolated','RSA Fixed','RSA Isolated'},'Location','best');
title('Superstructure Displacement')
xlabel('Displacement [in]');
ylabel('Height [ft]');


% Plot the Story shear forces
% Lets assume the isolator is 12" high and recalculate the drift
figure(4); hold on;
stories = 1:length(heights)-1;
plot(V_fixed, stories, '-o');
plot(V_isolated, [1 stories],'-o');
legend({'RSA Fixed','RSA Isolated','RSA Fixed','RSA Isolated'},'Location','best');
title('RHA Story Shear')
xlabel('Story shear [kips]');
ylabel('Story');


% Plot the Story shear forces
% Lets assume the isolator is 12" high and recalculate the drift
figure(5); hold on;
stories = 1:length(heights)-1;
plot(A_fixed, stories, '-o');
plot(A_isolated, [1 stories],'-o');
legend({'RSA Fixed','RSA Isolated','RSA Fixed','RSA Isolated'},'Location','best');
title('RHA Story Shear')
xlabel('Floor acceleration [g]');
ylabel('Story');

a=1






%% Question 5
% Initialization
nmodes = 5;
Gamma_n = NaN(nmodes,1);
u_n = NaN(9,nmodes,4251);
IDR = NaN(9,nmodes,4251);
f_n = NaN(9,nmodes,4251);
V_n = NaN(9,nmodes,4251);
phi = sphi;
phi_dr = [zeros(1,9); phi];


% Calculating for each mode
for i = 1:nmodes
    % Modal participation factor
    Gamma_n(i) = (phi(:,i)'*M*ones(size(phi(:,i))))/(phi(:,i)'*M*phi(:,i));
    % Displacement and acceleration histories for each mode
    [u,~,a,Sd,~,Sa,~,~] = RecurrenceSDOF(T(i),E,A,timestep,u0,v0,false);
    % Calculating for each time step
    for j = 1:4251
        % Displacement for each mode and timestep on all floors
        u_n(:,i,j) = Gamma_n(i)*phi(:,i).*u(j);
        % Calculating for each floor
        for k = 1:9
            % Interstory drift ratio for each floor for each timestep
            IDR(k,i,j) = (1/(Hi*12))*Gamma_n(i)*(phi_dr(k+1,i)-phi_dr(k,i)).*u(j);
        end
        % Lateral force for each mode, across all floors
        f_n(:,i,j) = Gamma_n(i)*phi(:,i)*mass*a(j)*g;
    end
    % Cumulatively summing the lateral forces in order to get shears
    for j = 1:4251
        fn_flip = flipud(f_n(:,i,j)); % Flipping because shear is the sum of forces from top to bottom
        V_n(:,i,j) = cumsum(fn_flip);
%         V_n(:,i,j) = flipud(V_n(:,i,j)); % Flipping again to get back into standard convention
    end
    for j = 1:4251
        V_n(:,i,j) = flipud(V_n(:,i,j));
    end
end





