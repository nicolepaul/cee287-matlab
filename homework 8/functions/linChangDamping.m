function [x] = linChangDamping(xo,T,damp)
% Adjust x with the damping modification factor proposed by Lin & Chang (2003)
    a = 1.303 + 0.436*log(damp); 
    B = 1 - a.*T.^0.3./(T+1).^0.65; 
    x = B.*xo'; 
end