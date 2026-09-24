clear
clc

A = [0 1; -24 -10];
[V,D] = eig(A);

% --- Input A ---
A = [0 1; -24 -10];

% --- Eigen-decomposition ---
[V,D] = eig(A);

disp('Eigenvalues (diag(D)):'), disp(diag(D))
disp('Eigenvectors (columns of V):'), disp(V)

% --- Verify A*V = V*D ---
residual = A*V - V*D;
fprintf('||A*V - V*D||_F = %.3e\n', norm(residual,'fro'));

% --- normalize eigenvectors to unit length ---
V = V ./ vecnorm(V);   

% If you want P and P^{-1} for diagonalization:
P  = V;
Pinv = inv(P);         
A_diag = Pinv*A*P;

function [V,D,residualNorm] = eig_with_check(A)
    [V,D] = eig(A);
    residualNorm = norm(A*V - V*D,'fro');
end