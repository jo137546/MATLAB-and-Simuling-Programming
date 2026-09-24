% ===== file: init_ss.m =====
% Define continuous-time plant
A = [0 1; -4 -3];
B = [0; 1];
C = [1 0];

% Make D the correct size (ny-by-nu) instead of scalar
D = zeros(size(C,1), size(B,2));

% Optional sample time (empty => continuous State-Space block)
Ts = [];   % e.g., Ts = 0.01 for discrete

% Robust sanity checks
nx = size(A,1);
    assert(ismatrix(A) && isnumeric(A) && size(A,2)==nx, 'A must be square (nx-by-nx).');
    assert(size(B,1)==nx, 'B must have nx rows.');
    assert(size(C,2)==nx, 'C must have nx columns.');
    assert(all(size(D) == [size(C,1) size(B,2)]), 'D must be ny-by-nu.');

% Optional deeper check
validate_ss(A,B,C,D);   % call the local function below

% ===== local function (at end of a script is OK in modern MATLAB) =====
function validate_ss(A,B,C,D)
    nx = size(A,1); ny = size(C,1); nu = size(B,2);
        assert(isequal(size(A),[nx nx]), 'A must be nx-by-nx.');
        assert(size(B,1)==nx, 'B must be nx-by-nu.');
        assert(size(C,2)==nx, 'C must be ny-by-nx.');
        assert(isequal(size(D),[ny nu]), 'D must be ny-by-nu.');
end