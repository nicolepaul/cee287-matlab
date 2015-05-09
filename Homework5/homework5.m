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
fprintf('----- Question 3 ------ \n')
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
% If we were using the code approach we would reduce by (R/Ie)
Cs = Sa;
nfloors = 9;
g = 386.1;
mass = 1800/g; % kips/g
stiffness = 1700; % kips/in
Hi = 118/9; %ft
[M, K] = computeMatrices(nfloors,mass,stiffness);
equivalentLateralForce(Tn,Cs,M,K,Hi)


%% Question 4
% Find the spectral acceleration of each mode
nmodes = 5;
Csm = zeros(nmodes,1);

% Calculate ehe spectral accleration, Csm for each mode
fprintf('Tm     Sa     Csm\n')
for i=1:5
    Tm = T(i);
    [~,~,~,~,~,Sa,~,~] = RecurrenceSDOF(Tm,E,A,timestep,u0,v0,false);
    Csm(i) = Sa;
    fprintf('%.3f   %.4f  %.4f  \n',Tm,Sa,Csm)
end

[F5,V5,U5,drift5] = Question4ModalAnalysis(nfloors,mass,stiffness,Csm,Hi);
[F3,V3,U3,drift3] = Question4ModalAnalysis(nfloors,mass,stiffness,Csm(1:3),Hi);

%close all;
figure; hold on;
floors = (9:-1:1)';

% Find the SS combination of first five modes
FnCombined = (sum(F5(:,1:5).^2,2));
VnCombined = (sum(V5(:,1:5).^2,2));
UxnCombined = (sum(U5(:,1:5).^2,2));
DxnCombined = (sum(drift5(:,1:5).^2,2));
clc
% Find the ratio of each mode contribution to total SS contribution
for mode=1:5
    fprintf('----- Mode %i -----',mode)
    FnContr = 100*(F5(:,mode).^2) ./ FnCombined
    VnContr = 100*(V5(:,mode).^2) ./ VnCombined
    UnContr = 100*(U5(:,mode).^2) ./ UxnCombined
    DnContr = 100*(drift5(:,mode).^2) ./ DxnCombined
end
    
    



for nmodes=1:5
    % Compute SRSS for 3 Modes
    FnCombined = (sum(F5(:,1:nmodes).^2,2)).^0.5;
    VnCombined = (sum(V5(:,1:nmodes).^2,2)).^0.5;
    UxnCombined = (sum(U5(:,1:nmodes).^2,2)).^0.5;
    DxnCombined = (sum(drift5(:,1:nmodes).^2,2)).^0.5;
    
    % Print forces along the height
    subplot(2,2,3); hold on;
    plot([FnCombined' 0], [floors' 0], '-o')
    title('Lateral Force')
    xlabel('Lateral Force [kips]')
    ylabel('Floor')
    grid on;
    axis([0 500 0 9]);
    
    % Print shear forces along the height
    subplot(2,2,4); hold on;
    plotSquare(flip(VnCombined), flip(floors),'-o');
    title('Shear Force')
    xlabel('Shear [kips]')
    ylabel('Floor')
    grid on;
    axis([0 2000 0 9]);
      
    % Print inelastic displacement along the height
    subplot(2,2,1); hold on;
    plot([UxnCombined' 0], [floors' 0],'-o')
    title('Lateral displacement, \delta_X')
    xlabel('\delta_X [in]')
    ylabel('Floor');
    grid on;
    axis([0 10 0 9]);
    
    % Print inelastic drift along the height
    subplot(2,2,2); hold on;
    plotSquare(flip(100*DxnCombined),flip(floors),'-o');
    title('Interstory Drift [%] \Delta_E')
    xlabel('\Delta_E [%]')
    ylabel('Floor');
    grid on;
    axis([0 0.6 0 9]); 
end
text = {'1 Mode','2 Modes','3 Modes','4 Modes','5 Modes'};
subplot(2,2,1); legend(text,'location','best');
subplot(2,2,2); legend(text,'location','best');
subplot(2,2,3); legend(text,'location','best');
subplot(2,2,4); legend(text,'location','best');



%% Question 5
nmodes = 5;
Gamma_n = NaN(nmodes,1);
u_n = NaN(9,nmodes,4251);
IDR = NaN(9,nmodes,4251);
f_n = NaN(9,nmodes,4251);
V_n = NaN(9,nmodes,4251);
phi = sphi;
phi_dr = [zeros(1,9); phi];
for i = 1:nmodes
    Gamma_n(i) = (phi(:,i)'*M*ones(size(phi(:,i))))/(phi(:,i)'*M*phi(:,i));
    [u,~,a,Sd,~,Sa,~,~] = RecurrenceSDOF(T(i),E,A,timestep,u0,v0,false);
    for j = 1:4251
        u_n(:,i,j) = Gamma_n(i)*phi(:,i).*u(j);
        for k = 1:9
            IDR(k,i,j) = (1/(Hi*12))*Gamma_n(i)*(phi_dr(k+1,i)-phi_dr(k,i)).*u(j);
        end
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
%%  Sum of modes for fl 9
floor=9;
figure;
subplot(5,1,1); plot(gmhw5(:,1),100*reshape(IDR(floor,1,:),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR [%]'); title(['Floor ' num2str(floor) ': Mode 1']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,2); plot(gmhw5(:,1),100*reshape(sum(IDR(floor,1:2,:),2),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 1 + 2']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,3); plot(gmhw5(:,1),100*reshape(sum(IDR(floor,1:3,:),2),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 1 + 2 + 3']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,4); plot(gmhw5(:,1),100*reshape(sum(IDR(floor,1:4,:),2),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 1 + 2 + 3 + 4']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,5); plot(gmhw5(:,1),100*reshape(sum(IDR(floor,1:5,:),2),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 1 + 2 + 3+ 4 + 5']); xlim([0 85]); ylim([-0.1 0.1]);

%% Each mode for fl 9
floor=9;
figure;
subplot(5,1,1); plot(gmhw5(:,1),100*reshape(IDR(floor,1,:),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR [%]'); title(['Floor ' num2str(floor) ': Mode 1']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,2); plot(gmhw5(:,1),100*reshape(IDR(floor,2,:),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 2']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,3); plot(gmhw5(:,1),100*reshape(IDR(floor,3,:),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 3']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,4); plot(gmhw5(:,1),100*reshape(IDR(floor,4,:),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 4']); xlim([0 85]); ylim([-0.1 0.1]);
subplot(5,1,5); plot(gmhw5(:,1),100*reshape(IDR(floor,5,:),1,4251)); grid on; xlabel('Time, t [s]'); ylabel('IDR  [%]'); title(['Floor ' num2str(floor) ': Mode 5']); xlim([0 85]); ylim([-0.1 0.1]);


%%
u_max = [0; max(abs(u_hist),[],2)];
IDR_max = [max(abs(IDR_hist),[],2); 0];
f_max = [0; max(abs(f_hist),[],2)];
V_max = [max(abs(V_hist),[],2); 0];


figure;
subplot(2,2,1); plot(u_max,0:9,'o-'); grid on; xlabel('Displacement [in]'); ylabel('Floor'); title('Maximum Displacements'); ylim([0 9]); xlim([0 6]);
subplot(2,2,2); plot(100*IDR_max,0:9,'o-'); grid on; xlabel('Interstory Drift Ratio [%]'); ylabel('Floor'); title('Maximum Drift Ratios'); ylim([0 9]); xlim([0 0.6]);
subplot(2,2,3); plot(f_max,0:9,'ro-'); grid on; xlabel('Lateral Force [kips]'); ylabel('Floor'); title('Maximum Lateral Forces'); ylim([0 9]); xlim([0 60]);
subplot(2,2,4); plot(V_max,0:9,'ro-'); grid on; xlabel('Story Shear [kips]'); ylabel('Floor'); title('Maximum Story Shear'); ylim([0 9]); xlim([0 400]);


%%

u_m = [ zeros(1,5); max(abs(cumsum(u_n,2)),[],3)];
IDR_m = [max(abs(cumsum(IDR,2)),[],3); zeros(1,5)];
f_m = [ zeros(1,5); max(abs(cumsum(f_n,2)),[],3)];
V_m = [max(abs(cumsum(V_n,2)),[],3); zeros(1,5)];

figure;
subplot(2,2,1); plot(u_m,0:9,'o-'); grid on; xlabel('Displacement [in]'); ylabel('Floor'); title('Maximum Displacements'); ylim([0 9]); xlim([0 6]); legend('1 mode','2 modes','3 modes','4 modes','5 modes','Location','best');
subplot(2,2,2); plot(100*IDR_m,0:9,'o-'); grid on; xlabel('Interstory Drift Ratio [%]'); ylabel('Floor'); title('Maximum Drift Ratios'); ylim([0 9]); xlim([0 0.6]); legend('1 mode','2 modes','3 modes','4 modes','5 modes','Location','best');
subplot(2,2,3); plot(f_m,0:9,'o-'); grid on; xlabel('Lateral Force [kips]'); ylabel('Floor'); title('Maximum Lateral Forces'); ylim([0 9]); xlim([0 60]); legend('1 mode','2 modes','3 modes','4 modes','5 modes','Location','best');
subplot(2,2,4); plot(V_m,0:9,'o-'); grid on; xlabel('Story Shear [kips]'); ylabel('Floor'); title('Maximum Story Shear'); ylim([0 9]); xlim([0 400]); legend('1 mode','2 modes','3 modes','4 modes','5 modes','Location','best');

