
clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lottery mean, covariance, and ellipses %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the six joint outcomes
% Column order: (1,A), (1,B), (1,C), (2,A), (2,B), (2,C)
%
% Row 1 contains the payout from scheme 1.
% Row 2 contains the payout from scheme 2.

payouts = [1,  1, 1, 2,  2, 2;
           3, -2, 3, 3, -2, 3];

% Probability of each corresponding column above
prob = [0.1; 0.2; 0.3; 0.2; 0.1; 0.1];

assert(abs(sum(prob) - 1) < 1e-12, ...
    'The outcome probabilities must sum to 1.');

%% Part 1: Expected payout vector
% Weighted sum of the six payout vectors
mu = payouts * prob;

fprintf('Expected payout for scheme 1: $%.2f\n', mu(1));
fprintf('Expected payout for scheme 2: $%.2f\n\n', mu(2));

%% Part 2: Second moment and centered covariance
M2 = zeros(2, 2);    % Uncentered second-moment matrix
P = zeros(2, 2);     % Covariance matrix about the mean

for k = 1:numel(prob)

    xk = payouts(:, k);    % Current two-component payout vector
    delta = xk - mu;       % Its deviation from the mean vector

    % Literal expression printed in the assignment
    M2 = M2 + prob(k) * (xk * xk.');

    % Covariance: use mean-centered payout vectors
    P = P + prob(k) * (delta * delta.');

end

disp('Uncentered second moment, M2:');
disp(M2);

disp('Centered covariance matrix, P:');
disp(P);

%% Part 3: Eigenvalues and eigenvectors
[V, D] = eig(P);

% Sort largest eigenvalue first, keeping eigenvectors matched
[lambda, order] = sort(diag(D), 'descend');
V = V(:, order);
D = diag(lambda);

disp('Eigenvalues, largest first:');
disp(lambda);

%% Generate and plot the covariance ellipses
theta = linspace(0, 2*pi, 500);

figure;
hold on;

for level = 1:3

    % Semi-axis lengths for this sigma level
    a = level * sqrt(lambda(1));
    b = level * sqrt(lambda(2));

    % Ellipse in the eigenvector coordinate frame
    eta = [a * cos(theta);
           b * sin(theta)];

    % Transform into the original payout coordinate frame
    xy = V * eta;

    % Translate the ellipse to the mean payout vector
    xy(1, :) = xy(1, :) + mu(1);
    xy(2, :) = xy(2, :) + mu(2);

    plot(xy(1, :), xy(2, :), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d-sigma', level));

end

% Mark the common center
plot(mu(1), mu(2), '+', ...
    'MarkerSize', 12, ...
    'LineWidth', 2, ...
    'DisplayName', 'Mean payout');

% Mark the actual payout pairs
% Some points overlap because A and C have the same payout.
plot(payouts(1, :), payouts(2, :), 'o', ...
    'MarkerSize', 7, ...
    'DisplayName', 'Possible payout pairs');

xlabel('Scheme 1 payout ($)');
ylabel('Scheme 2 payout ($)');
title('Lottery payout covariance ellipses');

axis equal;
grid on;
box on;
legend('Location', 'best');
hold off;