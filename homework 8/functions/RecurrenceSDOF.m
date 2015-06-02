function [u,v,a,Sd,Sv,Sa,PSv,PSa] = RecurrenceSDOF(Tn,E,A,dt,u0,v0,isPlot)
    %% 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RECURRENCE SDOF
    % Computes the seismic response of a linear elastic SDOF system using the
    % recurrence method
    % 
    % INPUTS:
    % Tn - undamped period of vibration of the system [s]
    % E - damping ratio of the system
    % A - ground motion acceleration time history [g]
    % u0 - initial displacement [in]
    % v0 - initial velocity [in/s]
    % isPlot - should we plot the response history
    % 
    % OUTPUTS:
    % u - relative displacement time history [cm]
    % v - relative velocity time history [in/s]
    % a - absolute acceleration time history [g]
    % Sd - spectral displacement [cm]
    % Sv - spectral velocity [in/s]
    % Sa - spectral acceleration [g]
    % PSv - pseudo velocity
    % PSa - pseudo acceleration
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Additional Inputs

    convg = 386;
    m=1; % This should cancel out, so set arbitrarily to 1


    %% Initial Computations

    % Determining natural and damped frequency [rad/s]
    wn = 2*pi/Tn;
    wd = wn*sqrt(1-E^2);
    k = wn^2*m;
    c = 2*E*m*wn;

    % Since dt is constant, calculation of constants below only needs to be
    % done once
    A1 = exp(-E*wn*dt)*(E*sin(wd*dt)/sqrt(1-E^2) + cos(wd*dt));
    B1 = exp(-E*wn*dt)*(sin(wd*dt)/wd);
    C1 = (1/k)*(2*E/(wn*dt) + exp(-E*wn*dt)*(((1-2*E^2)/(wd*dt) - E/sqrt(1-E^2))*sin(wd*dt) - (1 + 2*E/(wn*dt))*cos(wd*dt)));
    D1 = (1/k)*(1 - 2*E/(wn*dt) + exp(-E*wn*dt)*((2*E^2-1)*sin(wd*dt)/(wd*dt) + 2*E*cos(wd*dt)/(wn*dt)));
    A2 = -exp(-E*wn*dt)*(wn*sin(wd*dt)/sqrt(1-E^2));
    B2 = exp(-E*wn*dt)*(cos(wd*dt) - E*sin(wd*dt)/sqrt(1-E^2));
    C2 = (1/k)*(-1/dt + exp(-E*wn*dt)*((wn/sqrt(1-E^2) + E/(dt*sqrt(1-E^2)))*sin(wd*dt) + cos(wd*dt)/dt));
    D2 = (1/(k*dt))*(1 - exp(-E*wn*dt)*(E*sin(wd*dt)/sqrt(1-E^2) + cos(wd*dt)));

    %% Computing time histories

    % Determining force history
    p = -m*A*convg;
    n = numel(p);

    % Determining time steps
    t = (0:dt:n*dt-dt)';

    % Initialization of history variables
    u = zeros(n,1);
    v = zeros(n,1);
    a = zeros(n,1);

    % Including initial conditions
    u(1) = u0;
    v(1) = v0;
    a(1) = (p(1) - c*v0 - k*u0)/m;

    % Looping through all timesteps to calculate time history
    for i = 2:n % Initial conditions already given
        u(i) = A1*u(i-1) + B1*v(i-1) + C1*p(i-1) + D1*p(i);
        v(i) = A2*u(i-1) + B2*v(i-1) + C2*p(i-1) + D2*p(i);
        a(i) = (p(i) - c*v(i) - k*u(i))/m;
    end

    % Adding in ground acceleration and converting to g
    a = a/convg + A;

    % Determing spectral ordinates
    Sd = max(abs(u));
    Sv = max(abs(v));
    Sa = max(abs(a));

    % Determining pseudo spectral ordinates
    PSv = wn*Sd;
    PSa = wn^2*Sd/convg;

    % Plotting
    if isPlot
        figure; 
        subplot(4,1,1); plot(t,A,'k'); grid on; xlabel('Time [s]'); ylabel('Acceleration [g]'); title(['Time Histories for Tn = ' num2str(Tn) ', E = ' num2str(E)]); xlim([0 dt*n]);
        subplot(4,1,2); plot(t,u,'r'); grid on; xlabel('Time [s]'); ylabel('Displacement [cm]'); title('Displacement Time History');xlim([0 dt*n]);
        subplot(4,1,3); plot(t,v,'m'); grid on; xlabel('Time [s]'); ylabel('Velocity [cm/s]'); title('Velocity Time History'); xlim([0 dt*n]);
        subplot(4,1,4); plot(t,a,'b'); grid on; xlabel('Time [s]'); ylabel('Abs Acceleration [g]'); title('Absolute Acceleration Time History');xlim([0 dt*n]);
    end
end