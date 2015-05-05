%% Question 1a

% Compute the periods of vibration of the first three modes of vibration 
% (prepare  a  small  table  summarizing  periods,  frequencies  and
% circular  frequencies).
W = 1800; %kips/in
mass = W/386.4; % At each floor
stiffness = 1700; % At each floor
nfloors = 9;
[w,T,sphi,Gamma] = eigenvalueAnalysis(nfloors,mass,stiffness);
fprintf('Question 1:\n')
w
T
Gamma
printModalVectorTable(sphi);

figure; hold on; colorcell={'r','g','b'};
for i=1:3
    modeshape = sphi(:,i);
    gamma = Gamma(i);
    plot(gamma*[0 modeshape'],(0:nfloors),colorcell{i},'LineWidth',2);
end
grid on;
legend({'Mode 1','Mode 2','Mode 3'},'Location','best');
title('\Gamma_i\phi_i (First three modes)')
ylabel('Floor')
xlabel('\Gamma_i\phi_i')


%% Question 2

% Find the theoretical natural periods
fprintf('Question 2\n')
Tmodal = T/T(1);
Ttheory = 1./(2*(1:nfloors)'-1);
fprintf('Modal period ratios:\n')
disp(Tmodal)
fprintf('Theoretical period ratios:\n')
disp(Ttheory)

% Find the theoretical mode shapes
phiTheory = zeros(nfloors,3);
phiModal = zeros(nfloors,3);
for i = 1:3
    XoverH = (1:nfloors)'/nfloors;
    phiTheory(:,i) = (-1).^(i-1) * sin((2*i-1)*pi*XoverH/2);
    phiModal(:,i) = sphi(:,i)/sphi(end,i); % Roof normalize the modal results
end

fprintf('Actual Mode shapes:\n')
printModalVectorTable(phiModal)
fprintf('Theoretical mode shapes:\n')
printModalVectorTable(phiTheory)

% Find the theoretical modal participation factors
i = (1:3);
Gammat = 4*(-1).^(i-1)./(2*i*pi-pi);
fprintf('Actual Modal Participation Factors:\n');
disp(Gamma)
fprintf('Theoretical Modal Participation Factors:\n');
disp(Gammat)

% Plot the theoretical and real mode shapes
figure; hold on;
for i=1:3
    plot([0 phiModal(:,i)'], (0:nfloors),'g');
    plot([0 phiTheory(:,i)'],(0:nfloors),'b'); grid on;
end
legend({'Modal Analysis','Theoretical'});
title('Comparison of mode shapes (Roof normalized)')
ylabel('Floor')
xlabel('\phi_i (Roof normalized)')




%% Question 3
GM = dlmread('gm_hw_5.txt',' ',2,0);
timestep = GM(2,1) - GM(1,1);
conv_cm_to_g = 0.001019716;

% Spectra parameters
A = conv_cm_to_g*GM(:,2);
u0 = 0;
v0 = 0;
E = 0.035; % Damping ratio
periods = 0:0.02:3.0;
SA = zeros(length(periods),1);
SD = zeros(length(periods),1);

for i=1:length(periods)
   % Sd is returned in cm
   % Sa is returned in g
   Tn = periods(i);
   [~,~,~,Sd,~,Sa,~,~] = RecurrenceSDOF(Tn,E,A,timestep,u0,v0,false);
   SA(i) = Sa;
   SD(i) = Sd;
end

figure;
plot(periods,SA)
title('Linear elastic response spectra (\xi = 3.5%)')
xlabel('Period [s]')
ylabel('Spectral Acceleration [g]')

% Find the closest Sa ordinate to T1
% We could take this from the spectrum, or just recalculate
Tn = T(1);
[~,~,~,~,~,Sa,~,~] = RecurrenceSDOF(Tn,E,A,timestep,u0,v0,false);

% We use the Sa value from the spectra as our Cs
R = 8;
Ie = 1;
Cs = Sa/(R/Ie);
Cd = 5.5; % Special moment frame
nfloors = 9;
g = 386.1;
mass = 1800/g; % kips/g
stiffness = 1700; % kips/in
Hi = 118/9; %ft
[M, K] = computeMatrices(nfloors,mass,stiffness);
equivalentLateralForce(Tn,Cs,Cd,Ie,M,K,Hi)


%% Question 4
% Find the spectral acceleration of each mode
fprintf('Tm     Sa     Csm\n')
F = zeros(length(T),5);
Uxe = zeros(length(T),5);
Ux = zeros(length(T),5);
deltaX = zeros(length(T),5);
deltaXe = zeros(length(T),5);

for i=1:5
    Tm = T(i);
    [~,~,~,~,~,Sa,~,~] = RecurrenceSDOF(Tm,E,A,timestep,u0,v0,false);
    Csm = Sa/(R/Ie);
    An = Csm*g;
    fprintf('%.3f   %.4f  %.4f  \n',Tm,Sa,Csm)
    
    % Compute the modal forces for mode i
    F(:,i) = Gamma(i)*sphi(:,i)*mass*An;
    
    % Make floor 9 F(1,1)
    F = flipud(F);
    
    % Compute shear forces at each floor
    V(:,i) = cumsum(F(:,i));
    
    % Compute displacements at each floor (reduced)
    Uxe(:,i) = Gamma(i)*sphi(:,i)*An/(w(i))^2;
    Uxe(:,i) = flipud(Uxe(:,i));
   
    % Compute the amplifiied displacements
    Ux(:,i) = Cd*Uxe(:,i)/Ie;
    
    % Compute the interstorey drift (reduced)
    displacement = flipud(Uxe(:,i));
    deltaXe(:,i) = (displacement - [0 displacement(1:end-1)']')/(Hi*12);
    deltaXe(:,i) = flipud(deltaXe(:,i));
    
    % Compute the interstorey drift
    deltaX(:,i) = Cd*deltaXe(:,i)/Ie;
    
end

% Combine forces with
% SRSS-----------------------------------------------------------------------------------
FSRSS = (sum(F.^2, 2)).^0.5
FSRSS_with_floors = flipud(FSRSS)

%V
%Uxe
%Ux
deltaXe
deltaX






