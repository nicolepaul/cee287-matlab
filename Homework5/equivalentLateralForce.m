function equivalentLateralForce(Tn,Cs,Cd,I,M,K,Hi)
% Complete the equivalent lateral force procedure
    % @Tn - The natural period of the structure
    % @Cs - The seismic response coefficient
    % @Cd - ...
    % @I - The importance factor
    % @M - The mass matrix [kips/g]
    % @K - The stiffness matrix [kips/in]
    % @Hi - The height of each floor [ft]
    
    %Cs = SDS/(R/Ie)
    %Cs < SD1/(T*(R/Ie))       for T<Tl
    %Cs < SD1*TL/(T^2*(R/Ie))  for T>Tl
    %Cs > 0.044SDS*Ie >= 0.01
    
    % Calculate the weight of the structure
    % Assumes the mass matrix is in units kips/g
    % Based on Miranda, MDOF Seismic Analysis A, page 15
    W = sum(diag(M));
    V = Cs*W;
    fprintf('The seismic base shear, Vb = %.3f kips\n\n',V);
    
    
    heights = Hi*(1:length(M))';
    masses = diag(M);
    weights = masses; % Assumes mass in kip/g
    Fx = zeros(length(M),1);
    Vx = zeros(length(M),1);
    
    % Print the table of forces
    fprintf('Floor m[kip/g] wx[kip]     hx[ft]  Cvx   Fx[kip]    Vx[kip]\n')
    for i=1:length(M)
        floor = length(M)-i+1;
        k = getCVX_k(Tn);
        Cvx = weights(floor)*heights(floor)^k/(sum(weights.*heights.^k));
        Fx(floor) = Cvx*V;
        Vx(floor) = sum(Fx);
        fprintf('%i  %8.1f  %8.1f', floor, masses(floor), weights(floor));
        fprintf('%8.0f  %8.3f ', heights(floor), Cvx);
        fprintf('%8.3f  %8.3f \n', Fx(floor), Vx(floor));
    end
    
    deltaXE = K\Fx;
    deltaX = Cd/I*deltaXE;
    driftX = (deltaX - [0 deltaX(1:end-1)']')/Hi;
    
    % Print the table of displacements
    fprintf('\nFloor delXE[in]  delX[in]    IDR\n')
    for floor=1:length(M)
        floor = length(M)-floor+1;
        fprintf('%i  %10.3f ', floor, deltaXE(floor));
        fprintf('%10.3f %8.3f \n', deltaX(floor), driftX(floor));
    end
    
    % Print forces along the height
    figure;
    plot([0 Fx'], [0 heights'], '-o')
    title('Force Along Height')
    xlabel('Force [kips]')
    ylabel('Height [ft]')
    
    % Print shear forces along the height
    figure;
    plotSquare(Vx, heights,'-o');
    title('Shear Along Height')
    xlabel('Shear [kips]')
    ylabel('Height [ft]')
    
    % Print reduced elastic displacement along the height
    figure;
    plot([0 deltaXE'], [0 heights'],'-o')
    title('\delta_X_E Along Height')
    xlabel('\delta_X_E [in]')
    ylabel('Height[ft]')
    
    % Print inelastic displacement along the height
    figure;
    plot([0 deltaX'], [0 heights'],'-o')
    title('\delta_E Along Height')
    xlabel('\delta_E [in]')
    ylabel('Height[ft]')
    
    % Print inelastic drift along the height
    figure;
    plotSquare(driftX,heights,'-o');
    title('\Delta_E Along Height')
    xlabel('\Delta_E [in]')
    ylabel('Height[ft]')
end


function plotSquare(xvals,heights,options)
    % Same as line plot except everything is square
    % Assumes that xvals and heights are column vectors
    xvals = reshape(vertcat(xvals',xvals'),1,2*length(xvals));
    heights = reshape(vertcat(heights',heights'),1,2*length(heights));
    plot([xvals 0], [0 heights],options)
end


function k = getCVX_k(Tn)
    % Return the correct k for the calculation of CVX
    % k is a function of the period, T
    if Tn<=0.5
        k = 1;
    elseif Tn>=2.5
        k = 2;
    else
        k = 1+(Tn-0.5)/2;
    end
end