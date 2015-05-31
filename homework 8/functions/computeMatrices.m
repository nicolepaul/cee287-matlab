function [M, K] = computeMatrices(nfloors,mass,stiffness)
% Compute the stiffness and mass matrices for a nstorey structure
%   @nfloors: The number of storeys in the building
%   @stiffness: Either the stiffness at every floor (float) or a vector of stiffnesses
%   @mass: Either the mass at every floor (float) or a vector of masses
%   We use the convention that floor 1 is represented by K(1,1) ect
%   Where a vector is used for the mass or stiffness, the first item in the
%   vector should corrospond to the first floor.
    if length(stiffness)==1
        stiffness = stiffness*ones(1,nfloors);
    elseif length(stiffness)~=nfloors
        throw('Stiffness vector must be the same size as nfloors')
    end
    
    if length(mass)==1
        mass = mass*ones(1,nfloors);
    elseif length(mass)~=nfloors
        throw('Mass vector must be the same size as nfloors')
    end
    
    M = diag(mass);
    K = zeros(nfloors,nfloors);
    K(1,1) = stiffness(1);
        
    for i=2:nfloors
        kfloor = stiffness(i)*[1 -1; -1 1];
        K(i-1:i, i-1:i) = K(i-1:i, i-1:i) + kfloor;
    end
end

