function printModalVectorTable(sphi)
% Print a table of eigen vectors
%   By default our eigen vectors are organised so that sphi(1,1)
%   corrosponds to the first floor, first mode modeshape value.
%   Miranda arranges his tables so that the first floor is at the
%   bottom of the table. This function prints the matrix in this form
    fprintf('Floor ')    
    fprintf('  Mode %i  ',1:length(sphi))
    fprintf('\n')
    for row=1:length(sphi)
        floor = length(sphi) - row+1;
        fprintf('%4i ',floor)
        fprintf(' %8.3f ',sphi(floor,:))
        fprintf('\n')
    end
end

