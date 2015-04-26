function [u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAvgPrev(Tn, E, A, dt, u0, v0, uy, iflag)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AVERAGE ACCELERATION INTEGRATION
% Performs numerical integration using average acceleration approach
% Assumes acceleration input to be in g
% Calculates outputs in g, in/s, and in
% 
% Inputs:
% Tn - Undamped period of vibration [s]
% E - Damping ratio
% A - Ground motion acceleration time history [g]
% dt - Time step [s]
% u0 - Initial displacement [in]
% v0 - Initial velocity [in/s]
% uy - Yield displacement [in]
% iflag - Boolean set to 1 if iterations requested, 0 otherwise
% 
% Outputs:
% u - Relative displacement history [in]
% v - Relative velocity history [in/s]
% a - Absolute acceleration history [g]
% Sd - Spectral displacement ordinate [in]
% Sv - Spectral velocity ordinate [in/s]
% Sa - Spectral absolute acceleration ordinate [g]
% PSv - Pseudo velocity spectral ordinate [in/s]
% PSa - Pseudo acceleration spectral ordinate [g]
% Fs - Restoring force time history [kips]
% mu - Displacement ductility demand
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determining number of points given
M = 1; % assumed, should cancel out
g = 386; % assumed input is in [g]
P = -M*A*g;
N = numel(P);

% Setting up time vector
t = 0:dt:dt*N-dt;

% Determining necessary quantities
wn = 2*pi/Tn;
C = 2*E*wn*M;
K = wn^2*M;

% Initialization of output variables
u = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);
Fs = zeros(N,1);

% Initial conditions
u(1) = u0;
v(1) = v0;
a(1) = (P(1) - C*v0 - K*u0)/M;
uyn = -uy;
uyp = uy;      

% Performing numerical integration to get time history responses
for i = 2:N
    % Time histories
    u(i) = (P(i) + (M/dt^2)*(4*u(i-1)+4*v(i-1)*dt+a(i-1)*dt^2) + (C/dt)*(2*u(i-1)+v(i-1)*dt))/(4*M/dt^2 + 2*C/dt + K);
    v(i) = -v(i-1) + (2/dt)*(u(i)-u(i-1));
    a(i) = (4/dt^2)*(u(i)-u(i-1)-v(i-1)*dt-a(i-1)*dt^2/4);    
    % Inelastic check
    if u(i) > uyn && u(i) < uyp
       if v(i) <= 0
            Fs(i) = (K*uy + K*(u(i) - uyp));
       else
           Fs(i) = (-K*uy + K*(u(i) - uyn));
       end

    else
        if v(i) <= 0
            if v(i) < 0 && v(i-1) > 0
                uyp = u(i);
                uyn = u(i) - 2*uy;
            end
            Fs(i) = (K*uy + K*(u(i) - uyp));
            
        elseif v(i) > 0
            if v(i) > 0 && v(i-1) < 0
                uyn = u(i);
                uyp = u(i) + 2*uy;
            end
            Fs(i) = (-K*uy + K*(u(i) - uyn));
        end
        
        if u(i) < uyn
            Fs(i) = -K*uy;
        end
        
        if u(i) > uyp
            Fs(i) = K*uy;
        end
    
    end
    

end


% Get absolute acceleration
a = A + a/g;

% Determing spectral ordinates
Sd = max(abs(u));
Sv = max(abs(v));
Sa = max(abs(a));

% Determining pseudo spectral ordinates
PSv = wn*Sd;
PSa = wn^2*Sd/g;

% Calculating ductility demand
mu = Sd/uy;

% Creating figure to see results
set(0, 'defaultFigureColor', [1 1 1], 'defaultTextColor', [0 0 0]);
figure;
subplot(5,1,1);  plot(t,A,'k-'); grid on; xlabel('Time, t [s]'); ylabel('Acceleration [g]');  xlim([0 max(t)]); title(['Response History Plots, \mu = ' num2str(mu)],'FontName','Helvectica','FontWeight','bold','FontSize',14);
subplot(5,1,2); plot(t,u,'b-'); grid on; xlabel('Time, t [s]'); ylabel('Displacement [in]');  xlim([0 max(t)]); 
subplot(5,1,3); plot(t,v, 'r-'); grid on; xlabel('Time, t [s]'); ylabel('Velocity [in/s]'); xlim([0 max(t)]);
subplot(5,1,4); plot(t,a, 'g-'); grid on; xlabel('Time, t [s]'); ylabel('Acceleration [g]'); xlim([0 max(t)]);
subplot(5,1,5); plot(t,Fs, 'm-'); grid on; xlabel('Time, t [s]'); ylabel('Force [kips]'); xlim([0 max(t)]);

figure;
plot(u,Fs); grid on; xlabel('Displacement [in]'); ylabel('Resisting Force, F_s [kips]');
title('Force-Displacement');