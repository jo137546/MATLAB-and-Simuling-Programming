%%%%%%% Begin Boom Calculations %%%%%%%

% Define cell areas (if needed for later structural analysis)
CI = 1000 * (10^3);
CII = 260 * (10^3);
CIII = 180 * (10^3);

CI_prime = 40 * (10^3);
CII_prime = 105 * (10^3);
CIII_prime = 75 * (10^3);

% Boom areas
B = [600, 800, 800, 800, 800, 600]; % Effective areas of the booms (mm^2)

% Boom centroids (x-coordinates for the ROOT section)
x_boom_root = [590, 0, -790, -790, 0, 590]; % mm

% Boom centroids (x-coordinates for the FRONT section)
x_boom_front = [461.74, 0, -618.26, -618.26, 0, 461.74]; % mm (adjusted for taper)

% Boom centroids (y-coordinates, symmetric about the xz plane)
y_booms = [0, 0, 0, 0, 0, 0]; % mm (y-coordinates are zero due to symmetry)
 
% Root Section Centroid Calculation
disp('Root Section Centroid Calculation:');

% Initialize variables for weighted sums
weighted_sum_x_root = 0;
weighted_sum_y_root = 0;

% Total area of booms (same for both root and front)
total_area = sum(B);

% Loop through the booms to calculate weighted sums for the root section
for i = 1:length(B)
    weighted_sum_x_root = weighted_sum_x_root + B(i) * x_boom_root(i);
    weighted_sum_y_root = weighted_sum_y_root + B(i) * y_booms(i);
end

% Calculate root centroid coordinates
centroid_x_root = weighted_sum_x_root / total_area;
centroid_y_root = weighted_sum_y_root / total_area;

% Display results for the root section
disp(['Centroid (x) for Root: ', num2str(centroid_x_root), ' mm'])
disp(['Centroid (y) for Root: ', num2str(centroid_y_root), ' mm'])
 
 % Front Section Centroid Calculation
 disp('Front Section Centroid Calculation:');
 
 % Initialize variables for weighted sums
 weighted_sum_x_front = 0;
 weighted_sum_y_front = 0;
 
 % Loop through the booms to calculate weighted sums for the front section
 for i = 1:length(B)
    weighted_sum_x_front = weighted_sum_x_front + B(i) * x_boom_front(i);
    weighted_sum_y_front = weighted_sum_y_front + B(i) * y_booms(i);
end
 
% Calculate front centroid coordinates
centroid_x_front = weighted_sum_x_front / total_area;
centroid_y_front = weighted_sum_y_front / total_area;

% Display results for the front section
disp(['Centroid (x) for Front: ', num2str(centroid_x_front), ' mm'])
disp(['Centroid (y) for Front: ', num2str(centroid_y_front), ' mm'])