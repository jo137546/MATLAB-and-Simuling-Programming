%%%%%%%% Specify Wing Properties %%%%%%%

% Shear and Young's moduli (Pa)
G = 25 * (10^9); % Shear modulus
E = 70 * (10^9); % Young's modulus

% Thicknesses (m)
t1 = 2 / 1000;
t2 = 1.5 / 1000;
t3 = 1 / 1000;

% Lengths (m)
boom_length = 1000 / 1000;

l_12 = 600 / 1000;
l_23 = 800 / 1000;
l_34leading = 1200 / 1000;
l_45 = 800 / 1000;
l_56 = 600 / 1000;
l_16 = 210 / 1000;

l_34 = 320 / 1000;
l_25 = 320 / 1000;

l_12prime = 480 / 1000;
l_23prime = 620 / 1000;
l_34prime_leading = 750 / 1000;
l_45prime = 620 / 1000;
l_56prime = 480 / 1000;
l_16prime = 600 / 1000;

l_34prime = 200 / 1000;
l_25prime = 200 / 1000;

% Chord Lengths (m)
cl12 = 600 / 1000;
cl23 = 800 / 1000;
cl34leading = 750 / 1000;

cl12prime = 480 / 1000;
cl23prime = 620 / 1000;
cl34prime_leading = 750 / 1000;

%%%%%%centroid calculations%%%%%

% Areas (m^2)
A12 = l_12 * t3;
A23 = l_23 * t3;
A34 = l_34 * t1;
A45 = l_45 * t3;
A56 = l_56 * t3;
A16 = l_16 * t2;
A25 = l_25 * t2;
A34_leading = l_34leading * t3;

A12prime = l_12prime * t3;
A23prime = l_23prime * t3;
A34prime = l_34prime * t1;
A45prime = l_45prime * t3;
A56prime = l_56prime * t3;
A16prime = l_16prime * t2;
A25prime = l_25prime * t2;
A34prime_leading = l_34prime_leading * t3;

% x-Centroids
x_centroid = [
            590 / 2; % x12
    (790 / 2) + 590; % x23
          790 + 590; % x34
    (790 / 2) + 590; % x45
            590 / 2; % x56
                  0; % x16
                590; % x25
    sqrt(max(0, (2 * (cl34leading / pi)^2) - (320 / 2)^2)); % x34leading
];

x_centroid_prime = [
            470 / 2; % x12
    (610 / 2) + 470; % x23
          610 + 470; % x34
    (610 / 2) + 470; % x45
            470 / 2; % x56
                  0; % x16
                470; % x25
    sqrt(max(0, (2 * (cl34prime_leading / pi)^2) - (210 / 2)^2)); % x34leading
];

% y-centroids
y_centroid = 0; % symmetric
y_centroid_prime = 0; % symmetric


%initialize data points for x centroids
areas = [A12, A23, A34, A45, A56, A16, A25, A34_leading]; % Areas
x_centroids = [590/2, (790/2) + 590, 790 + 590, (790/2) + 590, 590/2, 0, 590, sqrt(max(0, (2 * (cl34leading / pi)^2) - (320 / 2)^2))]; % x-coordinates

numerator = 0; % Initialize the numerator (sum of A_i * x_i)
denominator = 0; % Initialize the denominator (sum of A_i)

for i = 1:length(areas)
    numerator = numerator + areas(i) * x_centroids(i); % Add A_i * x_i to numerator
    denominator = denominator + areas(i); % Add A_i to denominator
end

x_centroid_result = numerator / denominator; % Final centroid calculation

% Display Result %

disp(['The computed centroid (x-coordinate) is: ', num2str(x_centroid_result), ' m']);

%initialize data points for x centroids prime
areas = [A12prime, A23prime, A34prime, A45prime, A56prime, A16prime, A25prime, A34prime_leading]; % Areas
x_centroids_prime = [470 / 2, (610 / 2) + 470, 610 + 470, (610 / 2) + 470, 470 / 2, 0, 470, sqrt(max(0, (2 * (cl34prime_leading / pi)^2) - (210 / 2)^2))];

numerator = 0; % Initialize the numerator (sum of A_i * x_i)
denominator = 0; % Initialize the denominator (sum of A_i)

for i = 1:length(areas)
    numerator = numerator + areas(i) * x_centroids_prime(i); % Add A_i * x_i to numerator
    denominator = denominator + areas(i); % Add A_i to denominator
end

x_centroid_result = numerator / denominator; % Final centroid calculation

% Display Result %

disp(['The computed centroid (x-coordinate prime) is: ', num2str(x_centroid_result), ' m']);

% Display y centroid results

disp (['the computed centroid (y-coordinate) is:', num2str(y_centroid), 'm']); 
disp (['the computed centroid (y-coordinate prime) is:', num2str(y_centroid_prime), 'm']); 

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
x_boom = [-530, 0, 880, 880, 0, -530]

% Boom centroids (y-coordinates)
y_booms = [0,0,0,0,0,0]; %symmetric about xz plane

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