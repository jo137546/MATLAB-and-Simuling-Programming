clc; clear;

rng(1); % repeatable random noise

% Choose the TRUE polynomial coefficients
% a = [a1 a2 a3 a4 a5 a6 a7]^T

a_truth = rand(7,1); %returns m x n matrixof random numbers

% Generate 200 x measurements from -3 to 3
x = linspace(-3,3,200);

% y(t) function
y_truth = a_truth(1)*x.^6 +  a_truth(2)*x.^5 +  a_truth(3)*x.^4 + a_truth(4)*x.^3 + a_truth(5)*x.^2 + a_truth(6)*x + a_truth(7);
    
% x. esnures appropriate data management of the 200 element vector
   

%%%% Generate noisy measurements %%%%

variance = 0.1;
sigma = sqrt(variance);

noise = sigma * randn(size(x));

y_observed = y_truth + noise;

%%%% Plot true function and noisy measurements %%%%

figure;

% Noisy discrete measurements
plot(x, y_observed, 'ko','MarkerSize', 3);

hold on;

% True polynomial
plot(x, y_truth, 'b-','LineWidth', 2);

xlabel('x');
ylabel('y');

title('True Polynomial and Noisy Measurements');

legend('Noisy Measurements','True Function');

grid on;
hold off;

%%%% Batch Least Squares %%%%

% Create H matrix
H = zeros(200,7);

for i = 1:length(x)

    xi = x(i);

    H(i,:) = [xi^6, xi^5, xi^4, xi^3, xi^2, xi, 1];

end

%size(H)

%convert measurements into column vector
y_squiggly = y_observed(:);

%size(y_squiggly)

%%%% Estimate coefficients %%%%

a_hat = ((H' * H) ^ (-1)) * H' * y_squiggly;

%a_truth 
%a_hat

%%%% Coefficient estimation error %%%%

a_error = a_hat - a_truth;

%a_error

%%%% Estimated polynomial from Batch Least Squares %%%%

y_hat = H * a_hat;

%size(y_hat)

%%%% Batch Least Squares Quality-of-Fit Plot %%%%

figure;

% Noisy discrete measurements
plot(x, y_observed, 'ko', 'MarkerSize', 3, 'DisplayName', 'Noisy Measurements');

hold on;

% Estimated best-fit polynomial
plot(x, y_truth, 'b--', 'LineWidth', 2, 'DisplayName', 'TruePolynomial');

plot(x, y_hat, 'r-','LineWidth', 2, 'DisplayName', 'Batch LS Estimate');

xlabel('x');
ylabel('y');

title('True vs Batch Least-Squares Polynomial Fit');

legend('Location','best');
grid on;
hold off;

%%%% Compare true and estimated coefficients %%%%

coefficient = (1:7)';

comparison = table(coefficient,a_truth, a_hat, a_error, 'VariableNames', {'Coefficient','True_Value','Estimated_Value','Error'});

disp(comparison);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sequyential Least Squares alpha and beta approach %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 0.1; % measuerement noise variance

W = 1/R; % scalar measured weight as in problem 1

alpha = 1e3; %Large Number
beta = zeros(7,1); %small 7x1 "vector"

I = eye(7);

%%%% first measurement, first initialization %%%%

H1 = H(1,:); %1x7
y1 = y_squiggly(1); %scalar

P_seq = ((1/alpha^2 )* I + H1' * W * H1) \ I;

a_seq = P_seq * ((1/alpha)*beta + H1'*W*y1);

%size(a_seq)
%size(P_seq)

%%%% Sequential Least Squares Updates %%%%

a_seq_history = zeros(7,length(x));
a_seq_history(:,1) = a_seq;

for k = 2:length(x)

    % Current measurement model row
    Hk = H(k,:);

    % Current noisy measurement
    yk = y_squiggly(k);

    % Sequential least-squares gain
    K = P_seq * Hk' / (Hk * P_seq * Hk' + R);

    % Measurement residual / innovation
    residual = yk - Hk * a_seq;

    % Update coefficient estimate
    a_seq = a_seq + K * residual;

    % Update covariance
    P_seq = (I - K * Hk) * P_seq;

    % Save coefficient estimate
    a_seq_history(:,k) = a_seq;

end

%%%% Compare results %%%%
;
a_truth
a_hat      % batch estimate
a_seq    % sequential estimate

batch_error = a_hat - a_truth;
seq_error   = a_seq - a_truth;

batch_vs_seq = a_seq - a_hat

coefficient = (1:7)';

%%%% comparein a table %%%%
comparison = table(coefficient, a_truth, a_hat,  a_seq,batch_error, seq_error, 'VariableNames', {'Coefficient','Truth','Batch','Sequential', ...
    'Batch_Error','Sequential_Error'});

disp(comparison);

%%%% polynomial variables %%%%
y_seq = H * a_seq;

%%%% Batch vs Sequential Fit %%%%

figure;

plot(x, y_observed, 'ko', ...
    'MarkerSize', 3, ...
    'DisplayName', 'Noisy Measurements');

hold on;

plot(x, y_hat, 'r-', ...
    'LineWidth', 2, ...
    'DisplayName', 'Batch LS');

plot(x, y_seq, 'b--', ...
    'LineWidth', 2, ...
    'DisplayName', 'Sequential LS');

xlabel('x');
ylabel('y');

title('Batch vs Sequential Least-Squares Fit');

legend('Location','best');
grid on;
hold off;