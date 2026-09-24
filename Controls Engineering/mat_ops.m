
function C = mat_ops(A, op, B, varargin)
%MAT_OPS  Multiply or divide matrices (numeric & symbolic).
%   C = MAT_OPS(A, '*', B)     -> matrix multiply (A * B)
%   C = MAT_OPS(A, '\', B)     -> left division: solve A*X = B  (mldivide)
%   C = MAT_OPS(A, '/',  B)    -> right division: solve X*B = A (mrdivide)

%     'Steps' (false)   If true and op == '*', prints entry-wise dot products.
%
%   Notes:
%   - Works with numeric and symbolic inputs (requires Symbolic Math Toolbox for sym).
%   - Auto-promotes to symbolic when either A or B is symbolic, to keep exact arithmetic.
%   - Avoids inv(); for rectangular or rank-deficient cases:
%       * numeric: uses A\B or A/B (least-squares / min-norm as per MATLAB).
%       * symbolic: falls back to pinv(A) or pinv(B) if needed.

    p = inputParser;
    addParameter(p, 'Steps', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    showSteps = p.Results.Steps;

    % Validate types
    mustBeNumericOrSym(A);
    mustBeNumericOrSym(B);

    % Normalize op
    if ~(ischar(op) || (isstring(op) && isscalar(op)))
        error('Operation must be one of "*", "/", "\".');
    end
    op = char(op);

    % Auto-promote to symbolic if either input is symbolic
    isSym = isa(A,'sym') || isa(B,'sym');
    if isSym
        if ~isa(A,'sym'); A = sym(A); end
        if ~isa(B,'sym'); B = sym(B); end
    end

    switch op
        case '*'
            % Dimension check
            if size(A,2) ~= size(B,1)
                error('A*B dimension mismatch: size(A,2)=%d must equal size(B,1)=%d.', size(A,2), size(B,1));
            end
            C = A * B;

            if showSteps
                prettyPrintProductSteps(A,B);
            end

        case '\'
            % Solve A*X = B
            if size(A,1) ~= size(B,1)
                error('A\\B size mismatch: size(A,1)=%d must equal size(B,1)=%d.', size(A,1), size(B,1));
            end

            if isSym
                % Symbolic: prefer backslash; if it fails (e.g., singular/rectangular),
                % fallback to pseudoinverse for minimum-norm.
                try
                    C = A \ B;
                catch
                    C = pinv(A) * B;  % Moore-Penrose pseudoinverse (symbolic)
                end
            else
                % Numeric: MATLAB chooses robust factorization (QR/SVD/etc.)
                C = A \ B;
            end

        case '/'
            % Solve X*B = A
            if size(A,2) ~= size(B,2)
                error('A/B size mismatch: size(A,2)=%d must equal size(B,2)=%d.', size(A,2), size(B,2));
            end

            if isSym
                % Symbolic: prefer mrdivide; fallback to pseudoinverse if needed.
                try
                    C = A / B;
                catch
                    C = A * pinv(B);  % minimum-norm solution
                end
            else
                C = A / B;
            end

        otherwise
            error('Operation must be "*", "/", or "\".');
    end
end

% ---------- helpers ----------
function mustBeNumericOrSym(x)
    if ~(isnumeric(x) || isa(x,'sym'))
        error('Inputs must be numeric or symbolic (sym).');
    end
end

function prettyPrintProductSteps(A,B)
    % Prints C_ij = sum_k A_ik * B_kj
    m = size(A,1); n = size(B,2); kmax = size(A,2);
    fprintf('Entry-wise multiplication steps for C = A*B:\n');
    for i = 1:m
        for j = 1:n
            terms = strings(1,kmax);
            for k = 1:kmax
                terms(k) = sprintf('(%s)*(%s)', toChar(A(i,k)), toChar(B(k,j)));
            end
            fprintf('  C(%d,%d) = %s\n', i, j, strjoin(terms, ' + '));
        end
    end
end

function s = toChar(x)
    if isa(x,'sym')
        s = char(x);
    else
        s = num2str(x);
    end
end