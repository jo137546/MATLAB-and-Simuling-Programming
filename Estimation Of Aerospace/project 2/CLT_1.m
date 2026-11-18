clc; clear; close all;

%%%%%%%%%%%%%%%%%%%
%%%% problem 1 %%%%
%%%%%%%%%%%%%%%%%%%

%%% variables %%%
n_iter = [5, 10, 20 ,50,100]; % number of random values per sum N
n_trials = 1000;
store_hist = 30; %for histogram

n_cases = numel(n_iter);% give me number of matrix size

%%% allocate storage for historgram %%%

% rows = trials; columns = different values of N
trial_sum = zeros(n_trials, n_cases);
column_names = cell(1, n_cases); %represents each value of N

for j = 1:n_cases % loop through different cases of N

    N = n_iter(j); %cases per function of j'th X
    column_names{j} = sprintf('N_%d', N);

    for trial = 1: n_trials %for each parse through case, generate a random trial number
        
        X = rand(1,N); % generate all random numbers
        trial_sum(trial, j) = sum(X); %store sum of this trial

    end

    %%% plot 1000 sums for this N %%%
    
    figure;
    histogram(trial_sum(:, j),store_hist);
    
    xlabel('Sum of N uniform random numbers');
    ylabel('Frequency (number of trials)');
    
    title(sprintf('N = %d: %d trials', N, n_trials));
    
    grid on
end

%%% Tabulate the stored sums %%%
Trial = (1:n_trials).';

results_table = array2table([Trial, trial_sum], ...
    'VariableNames', [{'Trial'}, column_names]);

disp('First 10 trials:');
disp(results_table(1:10, :));

%%%%%%%%%%%%%%%%%%%
%%%% problem 2 %%%%
%%%%%%%%%%%%%%%%%%%

n_points = 1000;   % Number of positions used to calculate each curve

% Store the horizontal positions and corresponding Gaussian densities
% Each column corresponds to one value of N
x_gaussian = zeros(n_points, n_cases);
pdf_gaussian = zeros(n_points, n_cases);

for j = 1:n_cases

    N = n_iter(j);           % Select N = 5, 10, 20, or 50

    % Theoretical parameters specified in the assignment
    mu = N / 2;
    sigma = sqrt(N / 12);

    % Create evenly spaced positions between 0 and N
    x = linspace(0, N, n_points).';

    % Evaluate the Gaussian probability density at those positions
    f = normpdf(x, mu, sigma);

    % Store the calculated values for this N
    x_gaussian(:, j) = x;
    pdf_gaussian(:, j) = f;

    % Print the parameters
    fprintf('N = %d: mu = %.4f, sigma = %.6f\n', N, mu, sigma);

    % Plot the theoretical Gaussian curve
    figure;
    plot(x, f, 'LineWidth', 2);

    xlabel('Sum of N uniform random numbers');
    ylabel('Probability density');
    title(sprintf('Gaussian approximation for N = %d', N));

    xlim([0 N]);
    grid on;

end