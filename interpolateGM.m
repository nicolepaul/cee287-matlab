function A = interpolateGM(A, dt_old, dt_new)
    % Reduce the timestep of a ground motion record using interpolation
    N = numel(A);
    t = 0:dt_old:dt_old*N-dt_old;
    tnew = 0:dt_new:dt_new*N-dt_new;
    A = transpose(interp1(t,A,tnew));
end