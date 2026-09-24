% Input Parameters
S_y = 14764.8975; % Shear force in N
I_xx = 14764.9e6; % Moment of inertia in mm^4 (provided or calculated)
t_webs = [2, 2, 1.5, 2, 1.5]; % Thicknesses of webs in mm
lengths = [320, 200, 210, 200, 320]; % Segment lengths in mm
y_dist = [160, 100, 105, -100, -160]; % Distance of webs from neutral axis in mm
A_cells = [1e5, 2.6e5, 1.8e5]; % Areas of cells in mm^2
G = 25e3; % Shear modulus in MPa

% Initialize shear flow matrix (basic shear flow q_b)
q_b = zeros(1, length(lengths));

% Calculate basic shear flow for each segment
for i = 1:length(lengths)
    q_b(i) = (S_y / I_xx) * y_dist(i) * t_webs(i) * lengths(i);
end

% Display basic shear flow
disp('Basic shear flow for each segment:');
disp(q_b);

% Boundary conditions for open cells
q_b(1) = 0; % q_b,12
q_b(2) = 0; % q_b,23

% Contribution of each cell to total shear flow
delta_I = sum(q_b(3:end)); % Cell I contributions
delta_II = sum(q_b([2, 4])); % Cell II contributions
delta_III = sum(q_b([1, 5])); % Cell III contributions

% Solve for qs0 for each closed cell
syms qs0_I qs0_II qs0_III

eq1 = (qs0_I * sum(lengths) - delta_I) / (2 * A_cells(1) * G) == 0;
eq2 = (qs0_II * sum(lengths) - delta_II) / (2 * A_cells(2) * G) == 0;
eq3 = (qs0_III * sum(lengths) - delta_III) / (2 * A_cells(3) * G) == 0;

% Solve the equations for qs0
qs0_sol = solve([eq1, eq2, eq3], [qs0_I, qs0_II, qs0_III]);

qs0_I = double(qs0_sol.qs0_I);
qs0_II = double(qs0_sol.qs0_II);
qs0_III = double(qs0_sol.qs0_III);

% Display results
disp('Constant shear flow contributions (qs0):');
disp(['qs0_I: ', num2str(qs0_I), ' N/mm']);
disp(['qs0_II: ', num2str(qs0_II), ' N/mm']);
disp(['qs0_III: ', num2str(qs0_III), ' N/mm']);

% Calculate total shear flow in each segment
q_total = q_b;
q_total(3) = q_b(3) + qs0_I;
q_total(4) = q_b(4) + qs0_II;
q_total(5) = q_b(5) + qs0_III;

% Display total shear flows
disp('Total shear flows in each segment (N/mm):');
disp(q_total);
