clear; clc;

%%%% be careful with these variables dummy! %%%%

t = 0:1:25; % SYNTAX -> var = initialize:steps:finalize

K_meas = [0.5, 0.25, 0.125, 0, 0]; %k squiggly
K_truth = [0.2, 0.1, 0.05, 0.0001, 0.0001]; % k hat

lambda_meas = [-0.15, -0.06, -0.03];
lambda_truth = [-0.1, -0.05, -0.025];

omega_meas = [0.26, 0.55, 0.95];
omega_truth = [0.25, 0.5, 1];

delta_meas = [0.01, 0.01, 0.01];
delta_truth = [0, 0, 0];

sigma = 0.0002;

%%%% Starting estimate components %%%% math math math ...
A1 = K_meas(1) * exp(lambda_meas(1)*t);
B1 = omega_meas(1)*t + delta_meas(1);

C1 = K_meas(2) * exp(lambda_meas(2)*t);
D1 = omega_meas(2)*t + delta_meas(2);

E1 = K_meas(3) * exp(lambda_meas(3)*t);
F1 = omega_meas(3)*t + delta_meas(3);

%%%% True components %%%%
A2 = K_truth(1) * exp(lambda_truth(1)*t);
B2 = omega_truth(1)*t + delta_truth(1);

C2 = K_truth(2) * exp(lambda_truth(2)*t);
D2 = omega_truth(2)*t + delta_truth(2);

E2 = K_truth(3) * exp(lambda_truth(3)*t);
F2 = omega_truth(3)*t + delta_truth(3);

%%%% Four functions %%%%
theta_meas = A1.*cos(B1) + C1.*cos(D1) + E1.*cos(F1) + K_meas(4);
psi_meas   = A1.*sin(B1) + C1.*sin(D1) + E1.*sin(F1) + K_meas(5);

theta_truth = A2.*cos(B2) + C2.*cos(D2) + E2.*cos(F2) + K_truth(4);
psi_truth   = A2.*sin(B2) + C2.*sin(D2) + E2.*sin(F2) + K_truth(5);

%%%% Measurements %%%%

rng(1); % repeatable random noise

theta_observed = theta_truth + sigma*randn(size(t)); %true values + noise
psi_observed   = psi_truth   + sigma*randn(size(t)); %true values + noise

    %%%% Plot measurements sanity check %%%
    figure;

    % Starting estimate
    plot(t, theta_meas, 'mo', 'MarkerSize', 6);

    hold on;   % Keep this curve on the figure

    % True response
    plot(t, theta_truth, 'bo', 'MarkerSize', 6);

        xlabel('Time (s)');
        ylabel('\theta (Pitch)');
    
        title('Pitch: Measured Estimate vs True Response');

    legend('Starting Estimate', 'Truth');

    grid on;
    hold off;

%%%% model parameters x %%%%

x_est = [K_meas, lambda_meas, omega_meas, delta_meas]';


%%%% (finally omg i was freaking out) Finding Weight (W) %%%%

R = sigma^2 * eye(52);          % covariance matrix
W = (1/sigma^2) * eye(52);      % weight matrix


%%%% Iteration setup %%%%

maxIter = 5;                    % maximum 5 corrections
tol = 1e-6;                     % reasonable relative convergence tolerance

J_history = zeros(maxIter+1,1); % save J0, J1, J2, ...
x_history = zeros(14,maxIter+1);

x_history(:,1) = x_est;


%%%% BEGIN ITERATIVE NONLINEAR LEAST SQUARES %%%%

for iter = 1:maxIter


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Residuals using CURRENT parameter estimate %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delta_theta = theta_observed - theta_meas;
    delta_psi   = psi_observed   - psi_meas;

    % turn into desired delta y error matrix

    delta_y = zeros(52,1);

    delta_y(1:2:end) = delta_theta(:); % pitch error odd rows
    delta_y(2:2:end) = delta_psi(:);   % yaw error even rows


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Jacobian matrix H %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    H = zeros(52,14);

    for i = 1:length(t)

        ti = t(i);

        row_theta = 2*i - 1;
        row_psi   = 2*i;

        for j = 1:3

            % Extract CURRENT parameters for mode j
            k      = x_est(j);
            lambda = x_est(5+j);
            omega  = x_est(8+j);
            delta  = x_est(11+j);

            % Useful repeated quantities
            exponential = exp(lambda*ti);
            cosine      = cos(omega*ti + delta);
            sine        = sin(omega*ti + delta);

            % Derivative with respect to k_j
            H(row_theta,j) = exponential*cosine;
            H(row_psi,j)   = exponential*sine;

            % Derivative with respect to lambda_j
            H(row_theta,5+j) = k*ti*exponential*cosine;
            H(row_psi,5+j)   = k*ti*exponential*sine;

            % Derivative with respect to omega_j
            H(row_theta,8+j) = -k*ti*exponential*sine;
            H(row_psi,8+j)   =  k*ti*exponential*cosine;

            % Derivative with respect to delta_j
            H(row_theta,11+j) = -k*exponential*sine;
            H(row_psi,11+j)   =  k*exponential*cosine;

        end

        % k4 and k5 derivatives
        H(row_theta,4) = 1;
        H(row_theta,5) = 0;

        H(row_psi,4) = 0;
        H(row_psi,5) = 1;

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Scalar Cost Function J %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J_old = delta_y' * W * delta_y;

    % Save initial J only on first pass
    if iter == 1
        J_history(1) = J_old;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% delta x least squares solution %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Same equation you had, but use MATLAB backslash
    delta_x = (H' * W * H) \ (H' * W * delta_y);

    x_new = x_est + delta_x;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Recalculate model using x_new %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k1 = x_new(1);
    k2 = x_new(2);
    k3 = x_new(3);
    k4 = x_new(4);
    k5 = x_new(5);

    lambda1 = x_new(6);
    lambda2 = x_new(7);
    lambda3 = x_new(8);

    omega1 = x_new(9);
    omega2 = x_new(10);
    omega3 = x_new(11);

    delta1 = x_new(12);
    delta2 = x_new(13);
    delta3 = x_new(14);


    theta_new = ...
        k1*exp(lambda1*t).*cos(omega1*t + delta1) + ...
        k2*exp(lambda2*t).*cos(omega2*t + delta2) + ...
        k3*exp(lambda3*t).*cos(omega3*t + delta3) + ...
        k4;

    psi_new = ...
        k1*exp(lambda1*t).*sin(omega1*t + delta1) + ...
        k2*exp(lambda2*t).*sin(omega2*t + delta2) + ...
        k3*exp(lambda3*t).*sin(omega3*t + delta3) + ...
        k5;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Compare calculations against new prediction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delta_theta_new = theta_observed - theta_new;
    delta_psi_new   = psi_observed   - psi_new;

    %%% rinse and repeat %%%

    delta_y_new = zeros(52,1);

    delta_y_new(1:2:end) = delta_theta_new(:);
    delta_y_new(2:2:end) = delta_psi_new(:);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Updated cost function %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    J_new = delta_y_new' * W * delta_y_new;


    %%%%%%%%%%%%%%%%%%%%%%
    %%%% Evaluate delta J
    %%%%%%%%%%%%%%%%%%%%%%

    delta_J = abs(J_new - J_old);

    % Relative change makes convergence easier to judge
    relative_delta_J = delta_J / max(abs(J_old),eps);


    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Save iteration data
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    J_history(iter+1) = J_new;
    x_history(:,iter+1) = x_new;


    % Display progress
    fprintf(['Iteration %d: J_old = %.6e, J_new = %.6e, ' ...
             'relative change = %.6e\n'], ...
             iter, J_old, J_new, relative_delta_J);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% NEW becomes CURRENT for next iteration %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_est = x_new;

    theta_meas = theta_new;
    psi_meas   = psi_new;


    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Convergence check
    %%%%%%%%%%%%%%%%%%%%%%%%%

    if relative_delta_J < tol
        fprintf('Converged after iteration %d\n',iter);
        break
    end

end



%%%%% END ITERATIVE NONLINEAR LEAST SQUARES %%%%%

%%%% Cost Function Convergence Plot %%%%

iterations_completed = iter;

figure;

semilogy(0:iterations_completed, ...
         J_history(1:iterations_completed+1), ...
         'o-', 'LineWidth', 2, 'MarkerSize', 7);

xlabel('Iteration');
ylabel('Cost Function J');
title('Nonlinear Least Squares Convergence');

grid on;

%%%% Final Pitch Quality-of-Fit Plot %%%%


%store final best-fit curves %%%%
theta_bestfit = theta_new;
psi_bestfit   = psi_new;

x_bestfit = x_new;

figure;

% Discrete noisy measurements
plot(t, theta_observed, 'ko', ...
    'MarkerSize', 6, ...
    'DisplayName', 'Observed Measurements');

hold on;

% Best-fit estimated curve
plot(t, theta_new, 'r-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Estimated Best Fit');

xlabel('Time (s)');
ylabel('\theta (Pitch)');
title('Pitch: Estimated Fit vs Observed Measurements');

legend('Location','best');
grid on;
hold off;

%%%% Final Yaw Quality-of-Fit Plot %%%%

figure;

plot(t, psi_observed, 'ko', ...
    'MarkerSize', 6, ...
    'DisplayName', 'Observed Measurements');

hold on;

plot(t, psi_new, 'r-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Estimated Best Fit');

xlabel('Time (s)');
ylabel('\psi (Yaw)');
title('Yaw: Estimated Fit vs Observed Measurements');

legend('Location','best');
grid on;
hold off;

%%%% superimposed plots and sample data %%%%

figure;

plot(t, theta_observed, 'ko', 'MarkerSize', 6);
hold on;

plot(t, theta_bestfit, 'r-', 'LineWidth', 2);

xlabel('Time (s)');
ylabel('\theta');
title('Pitch Best Fit vs Measurements');

legend('Discrete Measurements','Best-Fit Curve');
grid on;