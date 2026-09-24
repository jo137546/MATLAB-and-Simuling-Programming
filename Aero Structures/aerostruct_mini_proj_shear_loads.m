%%%%% Shear Flow Calculations for Multi-Cell Wing %%%%%

% Inputs
Sy = 14764.8975; % N (shear force in y-direction)
Sx = 0; % N (shear force in x-direction, assumed zero for now)
Ixx = 95.15e6 / 1000^4; % Moment of inertia in m^4
Iyy = 1416.28e6 / 1000^4; % Moment of inertia in m^4
Ixy = 0; % Product moment of inertia, assumed zero
G = 25e9; % Shear modulus in Pa

% Thicknesses (m)
t = [0.002, 0.0015, 0.0015, 0.001, 0.001, 0.002]; % Thickness for each segment

% Segment lengths (m)
ds = [0.480, 0.320, 0.320, 0.480, 0.750, 0.750]; % Lengths in meters

% Initialize variables
q_b = zeros(1, length(ds));
tx_integral = 0;
ty_integral = 0;

% Coordinates of the webs (x and y in meters)
x_coords = [0.590, 0, -0.790, -0.790, 0, 0.590]; % x-coordinates
y_coords = [0.105, 0.160, 0.160, -0.160, -0.160, -0.105]; % y-coordinates

% Compute q_b for each segment
for i = 1:length(ds)
    tx_integral = tx_integral + x_coords(i) * ds(i); % Approximate integral for tx
    ty_integral = ty_integral + y_coords(i) * ds(i); % Approximate integral for ty

    q_b(i) = -((Sy * Ixx - Sx * Ixy) / (Ixx * Iyy - Ixy^2)) * tx_integral ...
             - ((Sx * Iyy - Sy * Ixy) / (Ixx * Iyy - Ixy^2)) * ty_integral;
end

% Debug: Display q_b values
disp('q_b values (basic shear flow for each segment):');
disp(q_b);

% Compute q_s,0 for the closed section
numerator_qs0 = sum(q_b .* ds); % Integral of q_b * ds

% Adjust denominator to prevent division by zero
denominator_qs0 = sum(ds ./ (G * t)); % Integral of ds / (G * t)
qs0 = -numerator_qs0 / denominator_qs0; % Calculate q_s,0

% Compute the final shear flow q_s for each segment
q_s = q_b + qs0;

% Segment names
segments = {'q12', 'q23', 'q34', 'q45', 'q56', 'q61'};

% Display results for each segment
disp('Shear Flow in Multi-Cell Wing (N/m):');
for i = 1:length(q_s)
    fprintf('%s = %.2f N/m\n', segments{i}, q_s(i));
end

% Display q_s,0
disp(['q_s,0 = ', num2str(qs0), ' N/m']);

% Debugging Outputs
disp('Debugging Outputs:');
disp(['Numerator for q_s,0: ', num2str(numerator_qs0)]);
disp(['Denominator for q_s,0: ', num2str(denominator_qs0)]);
disp('Thickness (t) values:');
disp(t);
disp('Segment lengths (ds):');
disp(ds);
disp('Sum of ds / (G * t):');
disp(sum(ds ./ (G * t)));