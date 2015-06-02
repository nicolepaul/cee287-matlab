function [u,v,a,f,vs,Sd,Sv,Sa] = modalTimeHistoryAbsolute(M,T,Gamma,Phi,E,Acc,dt)
    % Run a modal time history analysis
    % Inputs
    %   @M     - The mass matrix ()
    %   @T     - The modal periods
    %   @Gamma - The participation factors
    %   @Phi   - The mode shapes 
    %   @E     - The damping ratio
    %   @Acc   - The acceleration time history (g)
    %   @timestep - Acceleration time history timestep
    %   @useRelative - if true we will use the relative accelertion Dddot(t)
    %                  to calculate the total acceleration demand
    % Outputs
    %   @u(:,j) the displacement time history at floor j
    %   @a(:,j) the acceleration time history at floor j
    %   @f(:,j) the force time history at floor j
    %   @vb(:,j) the base shear time history at floor j
    %   @Sd(:,j) the peak floor displacement at floor j
    %   @Sa(:,j) the peak floor acceleration at floor j
    
    u0 = 0;
    v0 = 0;
    gravity = 386.1;
    nfloors = length(M);
    
    %           time      mode        floor
    D = zeros(length(Acc), length(T));
    U = zeros(length(Acc), length(T), length(M));
    V = zeros(length(Acc), length(T), length(M));
    A = zeros(length(Acc), length(T), length(M));
    F = zeros(length(Acc), length(T), length(M));
    C = zeros(length(Acc), length(T), length(M));
    
    for i = 1:length(T)
        [ui,vi,ai,~,~,~,~,~] = RecurrenceSDOF(T(i),E(i),Acc,dt,u0,v0,false);
        D(:,i) = ui;
        % Recurance SDOF will give us ai(total) - total acceleration
        % Displacement/acceleration at each floor in the structure
        U(:,i,:) = kron(Gamma(i)*Phi(:,i)', ui);
        V(:,i,:) = kron(Gamma(i)*Phi(:,i)', vi);
        A(:,i,:) = kron(Gamma(i)*Phi(:,i)', ai);
        F(:,i,:) = kron(Gamma(i)*Phi(:,i)'*M, ai*gravity);
        
        % Calculate the correction factor sum(Gamma*phi*ug(t))
        % The equation takes the same form as the one for A(:,i,:)
        C(:,i,:) = kron(Gamma(i)*Phi(:,i)', Acc);
    end
    
    % Postprocessing
    % Calculate the shear force for all modes and times
    VS = cumsum(F,3,'reverse');
  
    u = squeeze(sum(U,2));
    v = squeeze(sum(V,2));
    a = squeeze(sum(A,2)) - squeeze(sum(C,2)) + kron(Acc,ones(1,nfloors));
    f = squeeze(sum(F,2));
    vs = squeeze(sum(VS,2));
    
    Sd = max(abs(u),[],1);
    Sv = max(abs(u),[],1);
    Sa = max(abs(a),[],1);
end