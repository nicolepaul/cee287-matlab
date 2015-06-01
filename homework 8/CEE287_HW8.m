%% CEE 287: Homework 8
% Max Ferguson, Nicole Paul

%% Part A

%syms alpha real positive
alpha =  0.056409983382519;

gravity = 386.1; % in/s^2
designSd1 = 1.0;
targetDrift = 0.5/100; % 0.5%

% Define the height of each floor in the building
% The vector is defined from the bottom to top, starting at 0'
nfloors = 9;
hroof = 118;
heights = linspace(0,hroof,nfloors+1);

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

a1 = 1.303 + 0.436*log(E1);
a2 = 1.303 + 0.436*log(E2);

B1 = 1 - a1*T1^0.3/(T1+1)^0.65;
B2 = 1 - a2*T2^0.3/(T2+1)^0.65;

An1 = B1*designSd1/T1;
An2 = B2*designSd1;

Uj1 = Gamma1*phi1*An1*gravity/w1squared;
Uj2 = Gamma2*phi2*An2*gravity/w2squared;

Beta1 = GammaStruct*sphi(end);
temp_sphi = [0; sphi];
Beta2 = max(hroof*(temp_sphi(2:end)-temp_sphi(1:end-1))./(heights(2:end)*sphi(end))');

u = sqrt(Uj1.^2+Uj2.^2);
Sd1 = abs(u(2) - u(1));
IDR_maxSDOF = Sd1/((2/3)*hroof*12);
IDR_maxMDOF = (2/3)*Beta1*Beta2*IDR_maxSDOF;

fprintf('The drift %.4f [%%]\n',100*IDR_maxMDOF)
%alpha_design = eval(solve(IDR_maxMDOF-targetDrift==0))

%% 

targetdrift = 0.005;
f = @(alpha) targetdrift - get_drift(alpha);
alpha_design = fsolve2(f, 0.05);

m = 1000;
alpha_range = linspace(0,2,m);
drift_range = NaN(size(alpha_range));
for i = 1:m
    drift_range(i) = get_drift(alpha_range(i));
end

figure;
plot(alpha_range, drift_range, 'b-', alpha_design, targetdrift, 'ro'); grid on;
xlabel('Alpha, \alpha'); ylabel('Max Interstory Drift Ratio');
title('Influence of Alpha on IDR_{max}');


%% Part B

% Scale up stiffness from 2DOF system
MDOF_mass = nfloors*mass;
MDOF_stiffness = ks * MDOF_mass/ms;
stiffness_fixed = stiffness*ones(1,nfloors);
stiffness_isolated = [MDOF_stiffness,stiffness_fixed];

mass_fixed = mass*ones(1,nfloors);
mass_isolated = mass*ones(1,nfloors+1);


% 1) Get the stiffness and mass matrix of the equivalent structure
[M_fixed, K_fixed] = computeMatrices(nfloors, mass_fixed, stiffness_fixed);
[M_isolated, K_isolated] = computeMatrices(nfloors+1, mass_isolated, stiffness_isolated);

% Determining fundamental periods
nmodes = 2;
[~,T_fixed,sphi_fixed,Gamma_fixed] = eigenvalueAnalysis(nfloors,nmodes,mass_fixed,stiffness_fixed);
[~,T_isolated,sphi_isolated,Gamma_isolated] = eigenvalueAnalysis(nfloors+1,nmodes,mass_isolated,stiffness_isolated);

% Determining C_sm (seismic coefficient for each mode)
Sds = 1.5;
Sd1 = 1.0;
Csm_fixed = min([Sds*ones(2,1) Sd1./T_fixed],[],2);
Csm_isolated = min([Sds*ones(2,1) Sd1./T_isolated],[],2);

% Performing response spectrum analysis
[F_fixed,V_fixed,U_fixed,drift_fixed] = Question4ModalAnalysis(nfloors,mass_fixed,stiffness_fixed,Csm_fixed,heights(2:end));
[F_isolated,V_isolated,U_isolated,drift_isolated] = Question4ModalAnalysis(nfloors+1,mass_isolated,stiffness_isolated,Csm_isolated,heights);

%  2) Comparison plot of Story displacements for isolated vs fixed base structure