function [A, B, C, Q, R, yd] = CW_Matrix
% CW_Matrix  Circular-orbit relative-motion model and initial LQ weights.
% State order: [x; y; z; vx; vy; vz] in the leader's rotating LVLH frame.
% Axes: x radially outward, y along track, z along orbit angular momentum.
% Position is in m, velocity in m/s, and control U is acceleration in m/s^2.
% For the 10 kg follower, physical thrust is F = 10*U in N.

mu = 3.986e14;                % Earth's gravitational parameter, m^3/s^2
a = 6898e3;                  % Leader's circular-orbit radius, m
n = sqrt(mu / a^3);           % Mean motion, rad/s

% Kinematics occupy rows 1-3; CW acceleration equations occupy rows 4-6.
A = [0,       0,  0,    1,     0,   0;
     0,       0,  0,    0,     1,   0;
     0,       0,  0,    0,     0,   1;
     3*n^2,   0,  0,    0,     2*n, 0;
     0,       0,  0,   -2*n,   0,   0;
     0,       0, -n^2,  0,     0,   0];

% Each acceleration input acts directly on its corresponding velocity.
B = [zeros(3, 3);
     eye(3)];

% The measured outputs are the three relative positions.
C = [eye(3), zeros(3, 3)];

% Initial controller weights. Tune after adding estimation and saturation.
% Q weights position tracking; the sample forms Q_care = C'*Q*C.
% R penalizes acceleration commands. These are control-cost weights.
Q = eye(3);
R = 1e9 * eye(3);

% Keep the sample's six-entry reference convention for Q_care*yd.
% The desired measured position is C*yd = [50; 50; 0] m.
yd = [50; 50; 0; 0; 0; 0];
end