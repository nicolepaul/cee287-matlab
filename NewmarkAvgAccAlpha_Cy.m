function [u, v, a, Sd, Sv, Sa, PSv, PSa, Fs, mu] = NewmarkAvgAccAlpha_Cy(Tn, E, A, dt, u0, v0, Cy, alpha)
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
% Cy - seismic yield coefficient
% alpha - post elastic displacement
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
uy = Cy*M*g/K;
Fy = uy*K; 

% Initialization of output variables
u = zeros(N,1);
v = zeros(N,1);
a = zeros(N,1);
Fs = zeros(N,1);

% Newmark specific parameters
max_iterations = 10;
tol = 10^-3;
gamma = 1/2;
beta = 1/4;

% Intial conditions
u(1) = u0;
v(1) = v0;
[Fs(1),ktj] = getSpringForce(0,0,u0,K,alpha,Fy);
a(2) = (P(1) - C*v(1) - Fs(1))/M;

A1 = 1/(beta*dt^2)*M + gamma/(beta*dt)*C;
A2 = 1/(beta*dt)*M + (gamma/beta-1)*C;
A3 = (1/(2*beta)-1)*M + dt*(gamma/(2*beta)-1)*C;

for i=1:length(t)-1
    % These need to be nonlinear
    Fs(i+1) = Fs(i);
    Fsj = Fs(i+1);
    %ktj = getTangentStiffness(Fsj,K,Fy,alpha);
    %ktj = getSpringForce(Fsj,uprev,u,k,alpha,fy) 
%     [~,ktj] = getSpringForce(Fsj,u0,u(1),K,alpha,Fy)
    
    Pj = P(i+1) + A1*u(i) + A2*v(i) + A3*a(i);

    for j=1:max_iterations
        Rj = Pj - Fsj - A1*u(i+1);
        if abs(Rj)<tol
             fprintf('.')
            break
        elseif j==10
            fprintf('\nMax iterations exceeded\n')
        end
        ktj_hat = ktj+A1;
        u(i+1) = u(i+1)+(Rj/ktj_hat);
        % Update the spring force and tangent stiffness
        %Fsj = getSpringForce(Fs(i+1), u(i), u(i+1), K, Fy);
        [Fsj,ktj] = getSpringForce(Fs(i+1),u(i),u(i+1),K,alpha,Fy);
    end

    Fs(i+1) = Fsj;
    v(i+1) = gamma/(beta*dt)*(u(i+1)-u(i)) + (1-gamma/beta)*v(i) + dt*(1-gamma/(2*beta))*a(i);
    a(i+1) = 1/(beta*dt^2)*(u(i+1)-u(i)) - 1/(beta*dt)*v(i) - (1/(2*beta)-1)*a(i);
end
   
% Get absolute acceleration
a = a/g + A;

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
fprintf('\n')
set(0, 'defaultFigureColor', [1 1 1], 'defaultTextColor', [0 0 0]);
figure;
subplot(4,1,1); plot(t,A,'k-'); grid on; xlabel('Time, t [s]'); ylabel('Ground Acceleration [g]');  xlim([0 max(t)]); title(['Response History Plots, \mu = ' num2str(mu)],'FontName','Helvectica','FontWeight','bold','FontSize',14);
subplot(4,1,2); plot(t,u, 'r-'); grid on; xlabel('Time, t [s]'); ylabel('Displacement [in]'); xlim([0 max(t)]);
subplot(4,1,3); plot(t,v, 'g-'); grid on; xlabel('Time, t [s]'); ylabel('Velocity [in/s]'); xlim([0 max(t)]);
subplot(4,1,4); plot(t,Fs/(M), 'b-'); grid on; xlabel('Time, t [s]'); ylabel('Acceleration [in/s^2]'); xlim([0 max(t)]);


figure;
plot(u,Fs); grid on; xlabel('Displacement [in]'); ylabel('Resisting Force, F_s [kips]');
title('Force-Displacement');
end



function [fs,kt] = getSpringForceXXXX(fprev,uprev,u,k,alpha,fy)
    % Return the spring force, Fs,  at displacement u
    % The previous spring force and displacement are used to calculate the
    % new spring force. The spring force is limited to fy.
    
    fs = fprev - k*(uprev-u);
    
    if abs(fs)<=fy
        kt = k;
        fs = sign(fs)*fy;
    else
        kt = alpha*k;
    end
    
    % Subtract off the stiffness degradation
    fs = fs + u*k*alpha;
end




function [fs,kt] = getSpringForce(fprev,uprev,u,k,alpha,fy)
    % Return the spring force, Fs,  at displacement u
    % The previous spring force and displacement are used to calculate the
    % new spring force. The spring force is limited to fy.
    
    fs = fprev - k*(uprev-u);
    
    % Limit the maximum force to yield
    fs = sign(fs) * min(abs(fs),fy);
    
    % Return the tangent stiffness at displacement u
    % The tangent stiffness is either k, or alpha (post-yield)
    if abs(fs)>=fy
        kt = alpha*k;
    else
        kt = k;
    end
    
    % Subtract off the stiffness degradation
    fprintf('%.4f ')
    fs = fs + u*k*alpha;
    fprintf('%.4f ')
end
