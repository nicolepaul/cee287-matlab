% Test the computeMatrices function against the example from
% Chopra Dynamics of Structures Page 432

% From floor 1->n
mass = [1; 1; 0.5];
stiffness = [1; 7/9; 3/9];

Mexpected = diag(mass);
Kexpected = 1/9*[16  -7  0;
                 -7  10 -3;
                  0  -3  3];
              
[M,K] = computeMatrices(3,mass,stiffness);

fprintf('\nExpected Matrix (M):\n')
disp(Mexpected)
fprintf('\nActual Matrix (M)\n')
disp(M)

fprintf('\nExpected Matrix (K):\n')
disp(Kexpected)
fprintf('\nActual Matrix (K)\n')
disp(K)

