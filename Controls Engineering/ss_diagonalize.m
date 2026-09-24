function [Ad,Bd,Cd,Dd,T,Ti,info] = ss_diagonalize(A,B,C,D,real_modal)
% SS_DIAGONALIZE  Diagonalize (or real block-diagonalize) a state-space.
%
% Usage:
%   [Ad,Bd,Cd,Dd,T,Ti,info] = ss_diagonalize(A,B,C,D)          % real blocks
%   [Ad,Bd,Cd,Dd,T,Ti,info] = ss_diagonalize(A,B,C,D,false)    % allow complex
%
% Inputs:
%   A,B,C,D      State-space matrices (B,C,D optional; D defaults to zeros)
%   real_modal   (optional, default = true)  -> build a REAL modal form:
%                - real eigenvalues -> 1x1 blocks
%                - complex pairs    -> 2x2 real rotation/decay blocks
%                If false, returns a complex diagonalization.
%
% Outputs:
%   Ad,Bd,Cd,Dd  Transformed system in modal coordinates (x = T z)
%   T,Ti         Similarity matrices, Ti = inv(T) numerically via backslash
%   info         Struct with fields:
%                  .eigs        eigenvalues of A
%                  .blockDiag   block-diagonal Lambda used
%                  .residual    ||A*T - T*Lambda||_F
%                  .isRealForm  true if real 1x1/2x2 blocks

    if nargin < 2 || isempty(B), B = zeros(size(A,1),0); end
    if nargin < 3 || isempty(C), C = zeros(0,size(A,1)); end
    if nargin < 4 || isempty(D), D = zeros(size(C,1),size(B,2)); end
    if nargin < 5, real_modal = true; end

    n = size(A,1);
    [V,L] = eig(A);
    lam = diag(L);

    % check diagonalizability (geometric mult. == n)
    if rank(V) < n
        error('A is not diagonalizable (eigenvectors are not full rank).');
    end

    if real_modal && any(abs(imag(lam)) > 1e-12)
        % Build a REAL modal basis: 1x1 for real lambdas, 2x2 for complex pairs
        used = false(n,1);
        T  = [];           % real basis
        La = [];           % block-diagonal (real)
        for k = 1:n
            if used(k), continue; end
            if abs(imag(lam(k))) < 1e-12
                T  = [T, real(V(:,k))];
                La = blkdiag(La, real(lam(k)));
                used(k) = true;
            else
                j = find(~used & abs(lam - conj(lam(k))) < 1e-8, 1);
                if isempty(j)
                    error('Conjugate partner for complex eigenvalue not found.');
                end
                v = V(:,k);
                p = real(v); q = imag(v);
                a = real(lam(k)); b = imag(lam(k));
                T  = [T, p, q];
                La = blkdiag(La, [a -b; b a]);
                used([k j]) = true;
            end
        end
        Lambda = La;
        isRealForm = true;
    else
        T = V;
        Lambda = L;
        isRealForm = isreal(Lambda);
    end

    Ti = T \ eye(n);           % better-conditioned than inv(T)
    Ad = Ti*A*T;               % should equal Lambda up to numeric noise
    Bd = Ti*B;
    Cd = C*T;
    Dd = D;

    info.eigs      = diag(Lambda);
    info.blockDiag = Lambda;
    info.residual  = norm(A*T - T*Lambda,'fro');
    info.isRealForm = isRealForm;
end