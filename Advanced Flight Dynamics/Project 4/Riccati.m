function [R, B, P, g_aux] = Riccati(yd)
% Calculate the constant LQ design once per run, then apply the reference.
% Project4_Main clears this cache so edits to CW_Matrix take effect each run.
persistent cachedR cachedB cachedP referenceMap defaultReference

if isempty(cachedP)
    [A, B, C, Q, R, defaultReference] = CW_Matrix;
    Q_care = C' * Q * C;

    % Solve A'*P + P*A - P*B*(R\B')*P + Q_care = 0.
    P = care(A, B, Q_care, R);
    referenceMap = (P*B*(R\B') - A') \ Q_care;

    cachedR = R;
    cachedB = B;
    cachedP = P;
end

if nargin == 0
    yd = defaultReference;
end

R = cachedR;
B = cachedB;
P = cachedP;
g_aux = referenceMap * yd(:);
end
