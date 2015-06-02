function drift = get_drift(alpha)

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

ks = sphi(:,1)'*K*sphi(:,1);
ms = sphi(:,1)'*M*sphi(:,1);
ki = alpha*ks;
mi = beta*ms;

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

hroof = hroof*12;
H = (2/3)*hroof;
heights = heights*12;
Beta1 = GammaStruct*sphi(end);
temp_sphi = [0; sphi];
Beta2 = max(hroof*(temp_sphi(2:end)-temp_sphi(1:end-1))./(heights(2:end)*sphi(end))');
% hroof = hroof*12;
% H = (2/3)*hroof;
% Beta1 = Gamma1*phi1(1);
% Beta2 = hroof*(phi1(1)-phi1(2))/(H*phi1(1));

u = sqrt(Uj1.^2+Uj2.^2);
Sd1 = abs(u(2) - u(1));
% Sd1 = B1*An1*gravity/w1squared;
IDR_maxSDOF = Sd1/H;%Sd1/((2/3)*hroof*12);
IDR_maxMDOF = Beta1*Beta2*Sd1/hroof;%(2/3)*Beta1*Beta2*IDR_maxSDOF;

drift = IDR_maxMDOF;
end