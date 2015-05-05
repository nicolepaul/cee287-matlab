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
load gmhw5;
timestep = gmhw5(2,1) - gmhw5(1,1);
conv_cm_to_g = 0.001019716;

% Spectra parameters
A = conv_cm_to_g*gmhw5(:,2);
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

%% Question 5

nmodes = 5;
Gamma_n = NaN(nmodes,1);
u_n = NaN(9,nmodes,4251);
IDR = NaN(9,nmodes,4251);
f_n = NaN(9,nmodes,4251);
V_n = NaN(9,nmodes,4251);
phi = sphi;
phi_dr = phi';
for i = 1:nmodes
    Gamma_n(i) = (phi(:,i)'*M*ones(size(phi(:,i))))/(phi(:,i)'*M*phi(:,i));
    [u,~,a,Sd,~,Sa,~,~] = RecurrenceSDOF(T(i),E,A,timestep,u0,v0,false);
    conv_cm_to_in = 0.393701;
    for j = 1:4251
        u_n(:,i,j) = Gamma_n(i)*phi(:,i).*u(j)*conv_cm_to_in;
        IDR(:,i,j) = (1/(Hi*12))*Gamma_n(i)*(phi_dr(:,i+1)-phi_dr(:,i)).*u(j)*conv_cm_to_in;
        f_n(:,i,j) = Gamma_n(i)*phi(:,i)*a(j)*g;
    end
    for j = 1:4251
        fn_flip = flipud(f_n(:,i,j));
        V_n(:,i,j) = cumsum(fn_flip);
        V_n(:,i,j) = flipud(V_n(:,i,j));
    end
end

u_hist = reshape(sum(u_n,2),9,4251);
IDR_hist = reshape(sum(IDR,2),9,4251);
f_hist = reshape(sum(f_n,2),9,4251);
V_hist = reshape(sum(V_n,2),9,4251);


% figure;
% plot(gmhw5(:,1),u_hist); grid on;
% xlabel('Time, t [s]'); ylabel('Displacement, u_i [in]');
% legend('1','2','3','4','5','6','7','8','9','Location','best');
% title('Displacement History');
% 
% 
% figure;
% plot(gmhw5(:,1),IDR_hist); grid on;
% xlabel('Time, t [s]'); ylabel('IDR');
% legend('1','2','3','4','5','6','7','8','9','Location','best');
% title('Interstory Drift Ratio History');

u_max = max(abs(u_hist),[],2);
IDR_max = max(abs(IDR_hist),[],2);
f_max = max(abs(f_hist),[],2);
V_max = max(abs(V_hist),[],2);

figure;
subplot(2,2,1); plot(u_max,1:9,'o-'); grid on; xlabel('Displacement [in]'); ylabel('Floor'); title('Maximum Displacements');
subplot(2,2,2); plot(IDR_max,1:9,'o-'); grid on; xlabel('Interstory Drift Ratio'); ylabel('Floor'); title('Maximum Drift Ratios');
subplot(2,2,3); plot(f_max,1:9,'ro-'); grid on; xlabel('Lateral Force [kips]'); ylabel('Floor'); title('Maximum Lateral Forces');
subplot(2,2,4); plot(V_max,1:9,'ro-'); grid on; xlabel('Story Shear [kips]'); ylabel('Floor'); title('Maximum Story Shear');


