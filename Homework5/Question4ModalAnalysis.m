function [F,V,U,drift] = Question4ModalAnalysis(nfloors,mass,stiffness,Csm,Hi)
    % Find the spectral acceleration of each mode
    % Inputs:
    %   @nfloors: the number of floors in the building
    %   @mass: The mass at each floor kip.s^2/in
    %   @stiffness: The stiffness of each floor kip/in
    %   @Csm: The seismic coefficient
    %   @Cd: The displacement amplification factor
    %   @Ie: Importance modification factor
    %   @Hi: The hieght of each floor in ft
    % Outputs:
    %   @F: The force at each floor. F(1) is top floor
    %   @V: The shear force at each storey. V(1) is top floor
    %   @Uxe: The reduced elastic displacement at each floor. Uxe(1) is top floor
    %   @Ux: The inelastic displacement at each floor. Ux(1) is top floor
    %   @deltaXe: The reduced elastic drift at each storey. deltaXe(1) is top floor
    %   @deltaX: The inelastic drift at each storey. deltaX(1) is top floor
    
    [w,T,sphi,Gamma] = eigenvalueAnalysis(nfloors,mass,stiffness);
    
    nmodes = length(Csm);
    
    g = 386.1;
    F = zeros(length(T),nmodes);
    V = zeros(length(T),nmodes);
    U = zeros(length(T),nmodes);
    drift = zeros(length(T),nmodes);

    fprintf('Tm     Sa     Csm\n')
    for i=1:length(Csm)
        An = Csm(i)*g;
        fprintf('%.3f   %.4f  %.4f  \n',T(i),Csm(i))

        % Compute the modal forces for mode i
        F(:,i) = Gamma(i)*sphi(:,i)*mass*An;

        % Make floor 9 F(1,1)
        F(:,i) = flipud(F(:,i));

        % Compute shear forces at each floor
        V(:,i) = cumsum(F(:,i));

        % Compute displacements at each floor (reduced)
        U(:,i) = Gamma(i)*sphi(:,i)*An/(w(i))^2;
        U(:,i) = flipud(U(:,i));

        % Compute the interstorey drift (reduced)
        displacement = flipud(U(:,i));
        drift(:,i) = (displacement - [0 displacement(1:end-1)']')/(Hi*12);
        drift(:,i) = flipud(drift(:,i));
    end
end
    
    
    
    