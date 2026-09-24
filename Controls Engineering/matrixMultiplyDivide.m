A = [ 1 2 3; 4 5 6];
B = [ 1 2 ; 3 4];

function result = matrixMultiplyDivide(A, B, operation)
    % Define input matrices A and B
    A = [1 2 3; 4 5 6];
    B = [1 2; 3 4];
    
    % Define the operation to perform
    operation = 'multiply'; % or 'divide' depending on the desired operation
    
    % Call the plyDivide function
    result = plyDivide(A, B, operation);
    % matrixMultiplyDivide Multiplies or divides two matrices A and B.
    % Inputs:
    %   operation - 'multiply' or 'divide'
    % Output:
    %   result - Resulting matrix after the operation

    % Validate input matrices
    if ~isnumeric(A) || ~isnumeric(B)
        error('Both A and B must be numeric matrices.');
    end

    % Perform the specified operation
    switch lower(operation)
        case 'multiply'
            % Check if multiplication is possible
            if size(A, 2) ~= size(B, 1)
                error('Inner dimensions must match for multiplication.');
            end
            result = A * B;
        case 'divide'
            % Check if division is possible
            if size(B, 1) ~= size(B, 2)
                error('Matrix B must be square for division.');
            end
            % Calculate the inverse of B
            B_inv = inv(B);
            result = A * B_inv;
        otherwise
            error('Operation must be either "multiply" or "divide".');
    end
end