clear; clc;

%% 1. Generate the same measurements as the supplied Task 1 script
rng(6414,'twister');
pTrue = [0.05; 4; 0.2; -0.5; 10; pi/2];
s0True = [2; 0];
qTrue = [s0True; pTrue];       % [x(0); xdot(0); p1; ...; p6]
sigma = 0.1;
tm = (1:3000)'/10;
t = (0:0.02:300)';

% Keep the truth integration settings identical to Task 1.
optsTruth = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.03);
[xtrue_m,~,solTrue] = propagateModel(qTrue,tm,optsTruth);
noise = sigma*randn(numel(tm),1);
ymeas = xtrue_m + noise;       % Generate ONCE, outside the estimation loop
Ztrue = deval(solTrue,t).';
Strue = Ztrue(:,1:2);

%% 2. GLSDC settings and starting guesses
% The assignment requires 0.8 times the true initial states and parameters.
% Its hint allows a closer starting guess if this attempt does not converge.
% 0.8*xdot(0) is still zero; the initial velocity remains an estimated unknown.
guessScales = [0.8, 0.99];
maxIterations = 50;
costTolerance = 1e-9;
stepTolerance = 1e-7;
optsEstimate = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.03, ...
    'Events',@stopDivergingTrajectory);

attempts = cell(numel(guessScales),1);
nAttempts = 0;
selectedAttempt = [];
for a = 1:numel(guessScales)
    q0 = guessScales(a)*qTrue;
    fprintf('\nGLSDC starting guess: %.2f times truth\n',guessScales(a));
    attempts{a} = runGLSDC(q0,tm,ymeas,sigma,optsEstimate, ...
        maxIterations,costTolerance,stepTolerance);
    attempts{a}.scale = guessScales(a);
    nAttempts = a;
    fprintf('Status: %s\n',attempts{a}.status);
    if attempts{a}.converged
        selectedAttempt = a;
        break
    end
end
attempts = attempts(1:nAttempts);

% Do not present an unconverged attempt as a completed estimate.
if isempty(selectedAttempt)
    error(['Neither starting guess converged. Inspect the iteration output ' ...
        'and try a closer value in guessScales, as allowed by the assignment.']);
end

%% 3. Estimates, residuals, and local linearized covariance
fit = attempts{selectedAttempt};
qEstimate = fit.q;
[yfit,H,solFit] = propagateModel(qEstimate,tm,optsEstimate);
residual = ymeas-yfit;
Jfinal = 0.5*sum((residual/sigma).^2);
residualRMSE = sqrt(mean(residual.^2));

% R = sigma^2*I; W = R^(-1). Avoid forming a 3000-by-3000 matrix.
% Cov(q) is approximately (H'*W*H)^(-1), evaluated at the final estimate.
% This is a LOCAL linearized uncertainty, not a global nonlinear guarantee.
A = H/sigma;
columnScale = sqrt(sum(A.^2,1));
As = A./columnScale;
[~,Rqr] = qr(As,0);
D = diag(1./columnScale);
B = D*(Rqr\eye(8));
P = B*B.';
oneSigma = sqrt(diag(P));

names = {'x(0)'; 'xdot(0)'; 'p1'; 'p2'; 'p3'; 'p4'; 'p5'; 'p6'};
initialRequired = 0.8*qTrue;
initialUsed = guessScales(selectedAttempt)*qTrue;
estimationError = qEstimate-qTrue;
% The phase parameter is periodic; report its error on [-pi,pi].
estimationError(8) = atan2(sin(estimationError(8)),cos(estimationError(8)));
estimateTable = table(names,qTrue,initialRequired,initialUsed,qEstimate, ...
    estimationError,oneSigma,'VariableNames', ...
    {'Quantity','TrueValue','RequiredGuess','UsedGuess','Estimate','Error','OneSigma'});
disp(estimateTable);
fprintf('Selected starting scale: %.4f\n',guessScales(selectedAttempt));
fprintf('Accepted GLSDC updates: %d\n',numel(fit.Jhistory)-1);
fprintf('Final weighted cost J = %.10g\n',Jfinal);
fprintf('Residual RMS = %.10g; measurement standard deviation = %.3g\n', ...
    residualRMSE,sigma);
fprintf('Residual mean = %.10g\n',mean(residual));
fprintf('Reduced chi-square = %.10g\n',2*Jfinal/(numel(tm)-8));
fprintf('The p6 error is wrapped to [-pi,pi].\n');
fprintf(['Convergence means a local least-squares solution; it does not ' ...
    'guarantee the global minimum or exact recovery of every parameter.\n']);

%% 4. Ordinary cost plots for each attempted starting guess
figure('Name','Task 2 - GLSDC cost','NumberTitle','off','Visible','on','Color','w');
tiledlayout(1,nAttempts,'TileSpacing','compact');
for a = 1:nAttempts
    nexttile;
    J = attempts{a}.Jhistory;
    plot(0:numel(J)-1,J,'o-','LineWidth',1.2); grid on;
    xlabel('Accepted update number'); ylabel('Weighted least-squares cost, J');
    title(sprintf('Starting scale %.4f',attempts{a}.scale));
    if ~attempts{a}.converged
        text(0.03,0.95,attempts{a}.status,'Units','normalized', ...
            'VerticalAlignment','top','Interpreter','none');
    end
end

%% 5. Compare fitted motion with measurements and true motion
Zfit = deval(solFit,t).';
Sfit = Zfit(:,1:2);
figure('Name','Task 2 - Estimated motion','NumberTitle','off', ...
    'Visible','on','Color','w');
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(tm,ymeas,'.','Color',[0.7 0.7 0.7],'MarkerSize',4); hold on;
plot(t,Strue(:,1),'b','LineWidth',1);
plot(t,Sfit(:,1),'r--','LineWidth',1);
grid on; xlabel('Time, t'); ylabel('Position, x');
title('Measurements, true position, and estimated position');
legend('Measurements','True position','Estimated position');
nexttile;
plot(tm,ymeas,'.','Color',[0.6 0.6 0.6],'MarkerSize',6); hold on;
plot(t,Strue(:,1),'b','LineWidth',1.2);
plot(t,Sfit(:,1),'r--','LineWidth',1.2);
grid on; xlim([290 300]); xlabel('Time, t'); ylabel('Position, x');
title('Late motion, 290 to 300');
legend('Measurements','True position','Estimated position');

figure('Name','Task 2 - Measurement residuals','NumberTitle','off', ...
    'Visible','on','Color','w');
plot(tm,residual,'.','MarkerSize',5); hold on; yline(0,'k:');
grid on; xlabel('Time, t'); ylabel('Measured position - estimated position');
title('Final measurement residuals');
drawnow;

%% Local functions
function result = runGLSDC(q0,tm,y,sigma,opts,maxIter,costTol,stepTol)
    % Each iteration integrates the two states and all 16 sensitivities.
    result = struct('q',q0,'Jhistory',[],'qHistory',[], ...
        'converged',false,'status','Maximum iterations reached', ...
        'rejectedCost',NaN);
    q = q0;
    previousJ = Inf;
    for iteration = 0:maxIter
        try
            [ypred,H] = propagateModel(q,tm,opts);
        catch problem
            result.status = ['Integration failed: ' problem.message];
            return
        end
        e = y-ypred;                     % Measurements minus predicted position
        J = 0.5*sum((e/sigma).^2);
        fprintf('  iteration %2d: J = %.10g\n',iteration,J);

        % A full Gauss correction can fail when the starting guess is poor.
        % Stop this attempt and retain the preceding accepted point.
        if iteration>0 && J>previousJ*(1+1e-8)
            result.rejectedCost = J;
            result.status = 'Stopped because the next update increased the cost';
            return
        end
        result.q = q;
        result.Jhistory(end+1,1) = J;
        result.qHistory(:,end+1) = q;

        % Lecture update: dq=(H'*W*H)^(-1)*H'*W*e, W=I/sigma^2.
        % Backslash solves the equivalent weighted linear least-squares
        % problem without explicitly inverting the normal matrix.
        A = H/sigma;
        b = e/sigma;
        c = sqrt(sum(A.^2,1));
        if any(~isfinite(c)) || any(c==0)
            result.status = 'Stopped because a sensitivity column vanished';
            return
        end
        As = A./c;                       % Numerical column scaling only
        if rank(As)<8
            result.status = 'Stopped because the measurement matrix lost rank';
            return
        end
        scaledCorrection = As\b;
        dq = scaledCorrection./c.';
        relativeStep = max(abs(dq)./max(1,abs(q)));
        relativeCostChange = abs(previousJ-J)/max(previousJ,eps);
        if relativeStep<stepTol || (iteration>0 && relativeCostChange<costTol)
            result.converged = true;
            result.status = 'Converged';
            return
        end
        previousJ = J;
        q = q+dq;                        % Update all eight unknowns
    end
end

function [position,H,sol] = propagateModel(q,tm,opts)
    Phi0 = eye(2);
    Psi0 = zeros(2,6);
    z0 = [q(1:2); Phi0(:); Psi0(:)];
    sol = ode45(@(time,z) task2AugmentedODE(time,z,q(3:8)), ...
        [0 tm(end)],z0,opts);
    if sol.x(end)<tm(end)-1e-9
        error('The trial trajectory did not reach the last measurement time.');
    end
    Z = deval(sol,tm).';
    if any(~isfinite(Z),'all')
        error('The trial solution contains a nonfinite value.');
    end
    position = Z(:,1);
    % Column-major packing:
    % Phi(:) = [Phi11;Phi21;Phi12;Phi22].
    % Psi(:) = [Psi11;Psi21;Psi12;Psi22; ... ;Psi16;Psi26].
    % Only the POSITION row is measured.
    H = Z(:,[3 5 7 9 11 13 15 17]);      % 3000-by-8 measurement Jacobian
end

function dz = task2AugmentedODE(t,z,p)
    x = z(1);
    v = z(2);
    Phi = reshape(z(3:6),2,2);
    Psi = reshape(z(7:18),2,6);
    theta = p(5)*t+p(6);
    ds = [v; -p(1)*v-p(2)*x-p(3)*x^3-p(4)*sin(theta)];
    F = [0,1; -p(2)-3*p(3)*x^2,-p(1)];
    G = [zeros(1,6); ...
        -v,-x,-x^3,-sin(theta),-p(4)*t*cos(theta),-p(4)*cos(theta)];
    dPhi = F*Phi;
    dPsi = F*Psi+G;
    dz = [ds; dPhi(:); dPsi(:)];
end

function [value,isterminal,direction] = stopDivergingTrajectory(~,z)
    % Numerical safeguard only; no extra measurement or parameter prior.
    value = 1e3-abs(z(1));
    isterminal = 1;
    direction = -1;
end

