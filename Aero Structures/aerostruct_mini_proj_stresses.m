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