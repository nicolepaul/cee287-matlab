function equivalentLateralForce(Tn,Cs,M,K,Hi)
% Complete the equivalent lateral force procedure
    % @Tn - The natural period of the structure
    % @Cs - The seismic response coefficient
    % @M - The mass matrix [kips/g]
    % @K - The stiffness matrix [kips/in]
    % @Hi - The height of each floor [ft]
     
    % Calculate the weight of the structure
    % Assumes the mass matrix is in units kips/g
    % Based on Miranda, MDOF Seismic Analysis A, page 15
    
    g = 386.1; % in/s^2
    W = sum(diag(M))*g;
    V = Cs*W;
    fprintf('The seismic base shear, Vb = %.3f kips\n\n',V);
    
    floors = (1:length(M))';
    heights = Hi*(1:length(M))';
    masses = diag(M);
    weights = masses*g;
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
        fprintf('%8.0f  %9.3f ', heights(floor), Cvx);
        fprintf('%8.3f  %9.3f \n', Fx(floor), Vx(floor));
    end
    
    delta = K\Fx;
    drift = (delta - [0 delta(1:end-1)']')/(12*Hi);
    
    % Print the table of displacements
    fprintf('\nFloor displacement[in]   IDR[%%]\n')
    for i=1:length(M)
        floor = length(M)-i+1;
        fprintf('%i  %12.3f ', floor, delta(floor))
        fprintf('%12.3f\n',    100*drift(floor));
    end
    
    % Print inelastic displacement along the height
    figure;
    subplot(2,2,1);
    plot([0 delta'], [0 floors'],'-o')
    title('Lateral displacement, \delta_X')
    xlabel('\delta [in]')
    ylabel('Floor');
    grid on;
    axis([0 10 0 9]);
    
    % Print inelastic drift along the height
    subplot(2,2,2);
    plotSquare(100*drift,floors,'-o');
    title('Interstory Drift [%] \Delta/h')
    xlabel('\Delta/h [%]')
    ylabel('Floor');
    grid on;
    axis([0 2.0 0 9]);
    
    % Print forces along the height
    subplot(2,2,3);
    plot([0 Fx'], [0 floors'], '-o')
    title('Lateral Force')
    xlabel('Lateral Force [kips]')
    ylabel('Floor')
    grid on;
    axis([0 500 0 9]);
    
    % Print shear forces along the height
    subplot(2,2,4);
    plotSquare(Vx, floors,'-o');
    title('Shear Force')
    xlabel('Shear [kips]')
    ylabel('Floor')
    grid on;
    axis([0 2000 0 9]);
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