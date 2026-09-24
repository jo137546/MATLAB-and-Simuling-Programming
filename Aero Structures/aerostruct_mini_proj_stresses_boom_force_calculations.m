%%%%% Stress calculations %%%%%

% Calculated Moments of inertia
Mx = -2000; % Nm (bending moment around x-axis)
Ixx = 95.15e6 / 1000^4; % Convert mm^4 to m^4
Iyy = 1416.28e6 / 1000^4; % Convert mm^4 to m^4
Ixy = 0; % Product moment of inertia (assuming 0 for symmetry)

% Coordinates of booms relative to centroid (in mm)
x_booms = [590, 0, -790, -790, 0, 590]; % mm
y_booms = [105, 160, 160, -160, -160, -105]; % mm

% Convert coordinates to meters for calculation
x_booms_m = x_booms / 1000; % Convert to meters
y_booms_m = y_booms / 1000; % Convert to meters

% Calculate denominator
denominator = (Ixx * Iyy - Ixy^2);

% Calculate stress for each boom
sigma_z = zeros(1, length(x_booms)); % Preallocate stress array
for i = 1:length(x_booms_m)
    sigma_z(i) = Mx * (Iyy * y_booms_m(i) - Ixy * x_booms_m(i)) / denominator;
    disp(['Stress for boom ', num2str(i), ': ', num2str(sigma_z(i)), ' Pa']);
end

%%%%%% Perform boom angles / taper geometry %%%%%

% Deflections (in mm)
delx = [147.38, 27.38, -32.62, -32.62, 27.38, 147.38]; % mm
dely = [45, 60, 60, -60, -60, -45]; % mm

% Convert deflections to meters
delx_m = delx / 1000; % Convert to meters
dely_m = dely / 1000; % Convert to meters

% Deflection in z-direction (meters)
delz = 1000 / 1000; % Convert to meters

% Boom areas (in mm^2)
B = [600, 800, 800, 800, 800, 600]; % mm^2
B_m2 = B / 1e6; % Convert to m^2

% Initialize results matrix
results = zeros(length(x_booms), 4); % Columns: Sigma_z, Pz, Px, Py

% Force calculations
for i = 1:length(x_booms)
    % Calculate Pz
    Pz = B_m2(i) * sigma_z(i);

    % Calculate Px and Py
    Px = Pz * (delx_m(i) / delz);
    Py = Pz * (dely_m(i) / delz);

    % Store results
    results(i, :) = [sigma_z(i), Pz, Px, Py];
end

% Display results
disp('Boom Force Calculations:');
disp('Boom   Sigma_z (Pa)       Pz (N)      Px (N)      Py (N)');
for i = 1:length(x_booms)
    fprintf('%d      %.2e     %.2f     %.2f     %.2f\n', ...
        i, results(i, 1), results(i, 2), results(i, 3), results(i, 4));
end

%%%%% Reduced Shear Load Calculation %%%%%%

% Total shear force (N)
Sy = 14000;

% Extract Py from results
Py = results(:, 4); % Py is the 4th column of the results matrix

% Compute the reduced shear load
Sy_w = Sy - sum(Py);

% Display the reduced shear load
disp(['Reduced Shear Load (S_y,w): ', num2str(Sy_w), ' N']);

%%%%% General Shear Flow Calculation %%%%%

% Inputs
Sy = 14000; % N (shear force in y-direction)
Sx = 0; % N (shear force in x-direction, assumed zero for now)
Ixx = 95.15e6 / 1000^4; % Moment of inertia in m^4
Iyy = 1416.28e6 / 1000^4; % Moment of inertia in m^4
Ixy = 0; % Product moment of inertia, assumed zero

% Web thickness (assumed uniform for now)
t = 0.001; % Thickness in meters

% Coordinates of the webs (x and y in meters)
x_coords = [0.590, 0, -0.790, -0.790, 0, 0.590]; % in meters
y_coords = [0.105, 0.160, 0.160, -0.160, -0.160, -0.105]; % in meters

% Initialize variables for shear flow
qs = zeros(1, length(x_coords));

% Shear flow contributions
for i = 1:length(x_coords)
    % First integral: tx contribution
    tx_integral = 0;
    for j = 1:i
        tx_integral = tx_integral + x_coords(j) * t; % Approximate integral as sum
    end

    % Second integral: ty contribution
    ty_integral = 0;
    for j = 1:i
        ty_integral = ty_integral + y_coords(j) * t; % Approximate integral as sum
    end

    % Calculate shear flow using the formula
    qs(i) = -((Sy * Ixx - Sx * Ixy) / (Ixx * Iyy - Ixy^2)) * tx_integral ...
            - ((Sx * Iyy - Sy * Ixy) / (Ixx * Iyy - Ixy^2)) * ty_integral;
end

% Display the shear flow results
disp('Shear Flow in Each Web:');
for i = 1:length(qs)
    fprintf('Web %d: %.2f N/m\n', i, qs(i));
end