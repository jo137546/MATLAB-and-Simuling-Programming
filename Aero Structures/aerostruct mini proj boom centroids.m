%%%%%%begin boom calculations%%%%%%

%cell areas
CI = 1000*(10^3);
CII = 260*(10^3);
CIII = 180*(10^3);

CI_prime = 40*(10^3);
CII_prime = 105*(10^3);
CIII_prime = 75*(10^3);

% Boom areas
B = [600, 600, 800, 800, 800, 800];

% Boom centroids (x-coordinates)
x_boom = [530, 0, -880, -880,0, 530]

% Boom centroids (y-coordinates)
y_booms = [0,0,0,0,0,0] %symmetric about xz plane

%calculate total boom centroid:

% Initialize variables for weighted sums
weighted_sum_x = 0;
weighted_sum_y = 0;

% Total area of booms
total_area = sum(B);

% Loop to calculate the weighted sums
for i = 1:length(B)
    weighted_sum_x = weighted_sum_x + B(i) * x_boom(i);
    weighted_sum_y = weighted_sum_y + B(i) * y_booms(i);
end

% Calculate centroid coordinates
centroid_x = weighted_sum_x / total_area;
centroid_y = weighted_sum_y / total_area;

% Display results
disp(['Centroid (x): ', num2str(centroid_x), ' mm'])
disp(['Centroid (y): ', num2str(centroid_y), ' mm'])