function [w,T,sphi,Gamma] = eigenvalueAnalysis(nfloors,mass,stiffness)
% Compute periods of vibration, mode shapes and participation factors
%   @nfloors: The number of storeys in the building
%   @stiffness: Either the stiffness at every floor (float) or a vector of stiffnesses
%   @mass: Either the mass at every floor (float) or a vector of masses
%   We use the convention that floor 1 is represented by K(1,1) ect
%   Where a vector is used for the mass or stiffness, the first item in the
%   vector should corrospond to the first floor.
    [M,K] = computeMatrices(nfloors,mass,stiffness);
    [phi,w2] = eig(K,M);
    
    % Sort the modeshapes (Miranda, MDOF Structures Part A, page 34)
    [w2,index] = sort(diag(w2));
    for i=1:nfloors
        sphi(:,i) = phi(:,index(i));
    end
    
    % Obtaining the periods
    w = w2.^0.5;
    T = 2*pi./w;
   
    % Orthonormalizing the mode shapes
    sphi = roofNormalizeModeShapes(sphi,nfloors);
    %sphi = massNormalizeModeShapes(sphi,M,nfloors);
     
    % Calculate the modal participation factors
    Gamma = zeros(nfloors,1);
    for i=1:nfloors
        modeshape = sphi(:,i);
        Gamma(i) = modeshape'*M*ones(nfloors,1) / (modeshape'*M*modeshape);
    end
end


function sphi = massNormalizeModeShapes(sphi,M,nfloors)
    % Normalize the modeshapes such that sphi'*M*sphi = 1
    % See (Miranda, MDOF Structures Part A, page 35)
    for i=1:nfloors
       if sphi(:,i)'*M*ones(nfloors,1)>0
           sphi(:,i) = sphi(:,i)/sqrt((sphi(:,i)'*M*sphi(:,i)));
       else
            sphi(:,i) = -sphi(:,i)/sqrt((sphi(:,i)'*M*sphi(:,i)));
       end
    end
end


function sphi = roofNormalizeModeShapes(sphi,nfloors)
    % Normalize the modeshapes so that they have a roof displacement of 1
    % See (Miranda, MDOF Structures Part A, page 35)
    for i = 1:nfloors
        roof = sphi(end,i);
        sphi(:,i) = sphi(:,i)/roof;
    end
end

