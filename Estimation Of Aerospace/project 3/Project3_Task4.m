clear; clc;

Ntrials = 1000;               % Required number; reduce only for a test run
maxIterationsMC = 100;
initialScalesMC = [0.99 0.999 0.9999];
% 0.99 is the closer starting guess used by the successful reference fit.
% Retry the SAME measurement realization from a closer start if needed.
% Step halving is used when a full GLSDC correction increases the cost.
base = baselineTask2();
epochs = base.epochs;
[Plinear,~,~] = propagateErrorCovariance(base.qEstimate,base.P0, ...
    base.solEstimate,epochs,base.qTrue,base.trueEpoch);
rng(6415,'twister');          % Independent of the Task 3 reference realization
errors = nan(8,4,Ntrials);
qEstimatesMC = nan(8,Ntrials);
trialSuccess = false(Ntrials,1);
usedScale = nan(Ntrials,1);
trialCost = nan(Ntrials,1);
trialStatus = cell(Ntrials,1);
backtracked = false(Ntrials,1);
clockStart = tic;
for trial = 1:Ntrials
    ytrial = base.trueMeasurements + base.sigma*randn(numel(base.tm),1);
    fit = struct('converged',false,'status','No attempt completed');
    for a = 1:numel(initialScalesMC)
        fit = runDampedGLSDC(initialScalesMC(a)*base.qTrue,base.tm, ...
            ytrial,base.sigma,base.opts,maxIterationsMC,1e-9,1e-7);
        if fit.converged
            usedScale(trial) = initialScalesMC(a);
            break
        end
    end
    trialStatus{trial} = fit.status;
    if fit.converged
        trialSuccess(trial) = true;
        % Resolve equivalent amplitude-sign/frequency-sign/phase conventions
        % before comparing parameter errors across different fitted solutions.
        q = normalizeParameters(fit.q,base.qTrue(8));
        qEstimatesMC(:,trial) = q;
        ZE = deval(fit.sol,epochs);
        errors(1:2,:,trial) = ZE(1:2,:)-base.trueEpoch;
        pe = q(3:8)-base.qTrue(3:8);
        pe(6) = atan2(sin(pe(6)),cos(pe(6)));
        errors(3:8,:,trial) = repmat(pe,1,4);
        trialCost(trial) = fit.J;
        backtracked(trial) = fit.backtracked;
    end
    if mod(trial,10)==0 || trial==1 || trial==Ntrials
        elapsed = toc(clockStart);
        remaining = elapsed/trial*(Ntrials-trial);
        fprintf(['Trial %4d/%4d; converged %4d; elapsed %.1f min; ' ...
            'estimated remaining %.1f min\n'],trial,Ntrials, ...
            nnz(trialSuccess),elapsed/60,remaining/60);
        drawnow limitrate;
    end
end
Nconverged = nnz(trialSuccess);
failedTrialIndices = find(~trialSuccess);
if Nconverged<2
    error('Too few converged fits to compute a sample covariance.');
end
if Nconverged<Ntrials
    warning(['Only %d of %d realizations converged. Statistics below are ' ...
        'conditional on successful fits. Do not report 1000 converged ' ...
        'solutions unless Nconverged is 1000.'],Nconverged,Ntrials);
end
fprintf('\nConverged: %d of %d. Retries do not redraw measurement noise.\n', ...
    Nconverged,Ntrials);
fprintf('Successful fits requiring a closer restart: %d\n', ...
    nnz(usedScale>initialScalesMC(1)));
fprintf('Successful fits using at least one shortened correction: %d\n', ...
    nnz(backtracked));

%% Sample mean, centered covariance, and uncentered error second moment
meanError = zeros(8,4);
PmonteCarlo = zeros(8,8,4);
MSEmatrix = zeros(8,8,4);
sigmaLinear = zeros(8,4);
sigmaMC = zeros(8,4);
rmseMC = zeros(8,4);
coverage = zeros(4,3);
for k = 1:4
    E = reshape(errors(:,k,trialSuccess),8,Nconverged);
    meanError(:,k) = mean(E,2);
    Ec = E-meanError(:,k);
    PmonteCarlo(:,:,k) = (Ec*Ec.')/(Nconverged-1);
    MSEmatrix(:,:,k) = (E*E.')/Nconverged;
    sigmaLinear(:,k) = sqrt(diag(Plinear(:,:,k)));
    sigmaMC(:,k) = sqrt(diag(PmonteCarlo(:,:,k)));
    rmseMC(:,k) = sqrt(diag(MSEmatrix(:,:,k)));
    C = Plinear(1:2,1:2,k);
    distanceSquared = sum(E(1:2,:).*(C\E(1:2,:)),1);
    for level = 1:3
        coverage(k,level) = mean(distanceSquared<=level^2);
    end
end
covarianceDifference = PmonteCarlo-Plinear;
stdRatio = sigmaMC./sigmaLinear;
comparisonTable = table(kron(epochs(:),ones(8,1)),repmat(base.names,4,1), ...
    meanError(:),sigmaMC(:),sigmaLinear(:),stdRatio(:),rmseMC(:), ...
    'VariableNames',{'Time','Quantity','MeanError','MCStandardDeviation', ...
    'LinearStandardDeviation','MCOverLinearStd','MCRootMeanSquareError'});
disp(comparisonTable);
coverageTable = table(epochs(:),coverage(:,1),coverage(:,2),coverage(:,3), ...
    'VariableNames',{'Time','Inside1Sigma','Inside2Sigma','Inside3Sigma'});
disp(coverageTable);
fprintf('Joint Gaussian reference coverage: %.6f %.6f %.6f\n', ...
    1-exp(-[1 4 9]/2));
fprintf(['PmonteCarlo is centered on the sample mean. MSEmatrix estimates ' ...
    'E[error*error''] and includes bias.\n']);
for k = 1:4
    fprintf('\nt = %g: linear state covariance\n',epochs(k));
    disp(Plinear(1:2,1:2,k));
    fprintf('Monte Carlo centered state covariance\n');
    disp(PmonteCarlo(1:2,1:2,k));
    fprintf('State error second-moment matrix\n');
    disp(MSEmatrix(1:2,1:2,k));
    fprintf('Full 8-by-8 covariance difference, MC minus linear\n');
    disp(covarianceDifference(:,:,k));
    identity = (Nconverged-1)/Nconverged*PmonteCarlo(:,:,k) ...
        +meanError(:,k)*meanError(:,k).';
    assert(norm(MSEmatrix(:,:,k)-identity,'fro')< ...
        1e-10*max(1,norm(MSEmatrix(:,:,k),'fro')));
end

%% Native MATLAB plots
figure('Name','Task 4 - Monte Carlo error clouds','NumberTitle','off', ...
    'Visible','on','Color','w','Position',[80 80 1150 850]);
tiledlayout(2,2,'TileSpacing','compact');
for k = 1:4
    nexttile; hold on;
    E = reshape(errors(:,k,trialSuccess),8,Nconverged);
    hs = plot(E(1,:),E(2,:),'.','Color',[0.65 0.65 0.65],'MarkerSize',6);
    h = drawEllipses(Plinear(1:2,1:2,k),[0;0],[1 2 3],'-');
    hm = drawEllipses(PmonteCarlo(1:2,1:2,k),meanError(1:2,k),3,'--');
    hb = plot(meanError(1,k),meanError(2,k),'kx','LineWidth',1.5,'MarkerSize',9);
    grid on; axis equal; xlabel('Position error'); ylabel('Velocity error');
    title(sprintf('t = %g, %d converged trials',epochs(k),Nconverged));
    legend([hs;h(:);hm(:);hb],{'Monte Carlo errors','Linear 1 sigma', ...
        'Linear 2 sigma','Linear 3 sigma','MC 3 sigma','MC mean'}, ...
        'Location','best');
end

figure('Name','Task 4 - Standard deviation comparison','NumberTitle','off', ...
    'Visible','on','Color','w','Position',[90 90 1100 700]);
tiledlayout(2,2,'TileSpacing','compact');
for k = 1:4
    nexttile; bar(stdRatio(:,k)); hold on; yline(1,'k--');
    xticks(1:8); xticklabels(base.names); xtickangle(30);
    ylabel('Monte Carlo std / linear std'); grid on;
    title(sprintf('t = %g',epochs(k)));
end

figure('Name','Task 4 - Joint ellipse coverage','NumberTitle','off', ...
    'Visible','on','Color','w');
bar([coverage;1-exp(-[1 4 9]/2)]);
xticklabels({'t=0','t=100','t=200','t=300','Gaussian reference'});
ylabel('Fraction inside the joint ellipse'); ylim([0 1]); grid on;
legend('1 sigma','2 sigma','3 sigma','Location','best');
title(sprintf('Coverage from %d converged solutions',Nconverged));
drawnow;


%% Local functions

function base = baselineTask2()
    % Reproduce the single measurement realization and full-step Task 2 fit.
    base.qTrue = [2;0;0.05;4;0.2;-0.5;10;pi/2];
    base.sigma = 0.1;
    base.tm = (1:3000)'/10;
    base.epochs = [0 100 200 300];
    base.names = {'x';'xdot';'p1';'p2';'p3';'p4';'p5';'p6'};
    optsTruth = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.03);
    base.opts = odeset(optsTruth,'Events',@stopDivergingTrajectory);
    rng(6414,'twister');
    [base.trueMeasurements,~,solTrue] = propagateModel(base.qTrue,base.tm,optsTruth);
    base.y = base.trueMeasurements+base.sigma*randn(numel(base.tm),1);
    ZE = deval(solTrue,base.epochs);
    base.trueEpoch = ZE(1:2,:);
    for scale = [0.8 0.99]
        fprintf('\nReference Task 2 fit: starting scale %.4f\n',scale);
        fit = runGLSDC(scale*base.qTrue,base.tm,base.y,base.sigma, ...
            base.opts,50,1e-9,1e-7);
        fprintf('%s\n',fit.status);
        if fit.converged
            break
        end
    end
    if ~fit.converged
        error('The reference Task 2 fit failed. Inspect the iteration output.');
    end
    base.qEstimate = fit.q;
    [~,H,base.solEstimate] = propagateModel(fit.q,base.tm,base.opts);
    base.P0 = linearCovariance(H,base.sigma);
end

function P = linearCovariance(H,sigma)
    A = H/sigma;
    c = sqrt(sum(A.^2,1));
    if any(c==0) || any(~isfinite(c))
        error('Invalid measurement-sensitivity column.');
    end
    As = A./c;
    if rank(As)<8
        error('The local information matrix does not have rank eight.');
    end
    [~,R] = qr(As,0);
    B = diag(1./c)*(R\eye(8));
    P = B*B.';
    P = (P+P.')/2;
end

function [Pz,errors,T] = propagateErrorCovariance(q,P0,sol,epochs,qTrue,trueEpoch)
    % z(t) = [x(t); xdot(t); p1; ...; p6].
    % dz(t) = T(t)*dq(0), T = [Phi Psi; zeros(6,2) eye(6)].
    % The parameters are constant; their uncertainty still affects motion.
    ZE = deval(sol,epochs);
    Pz = zeros(8,8,numel(epochs));
    T = zeros(8,8,numel(epochs));
    errors = zeros(8,numel(epochs));
    for k = 1:numel(epochs)
        Phi = reshape(ZE(3:6,k),2,2);
        Psi = reshape(ZE(7:18,k),2,6);
        T(:,:,k) = [Phi,Psi;zeros(6,2),eye(6)];
        Pz(:,:,k) = T(:,:,k)*P0*T(:,:,k).';
        Pz(:,:,k) = (Pz(:,:,k)+Pz(:,:,k).')/2;
        errors(:,k) = [ZE(1:2,k)-trueEpoch(:,k);q(3:8)-qTrue(3:8)];
        errors(8,k) = atan2(sin(errors(8,k)),cos(errors(8,k)));
    end
end

function h = drawEllipses(C,center,levels,lineStyle)
    C = (C+C.')/2;
    [V,Lambda] = eig(C);
    eigenvalues = diag(Lambda);
    if any(eigenvalues<=0) || any(~isfinite(eigenvalues))
        error('The position-velocity covariance must be positive definite.');
    end
    angle = linspace(0,2*pi,301);
    unitCircle = [cos(angle);sin(angle)];
    L = V*diag(sqrt(eigenvalues));
    colors = lines(3);
    h = gobjects(numel(levels),1);
    for j = 1:numel(levels)
        radius = levels(j);
        points = center+radius*L*unitCircle;
        h(j) = plot(points(1,:),points(2,:),'LineStyle',lineStyle, ...
            'Color',colors(radius,:),'LineWidth',1.3);
    end
end

function q = normalizeParameters(q,phaseReference)
    % Equivalent sine representations give exactly the same physical forcing.
    if q(7)<0
        q(7) = -q(7);
        q(6) = -q(6);
        q(8) = -q(8);
    end
    if q(6)>0
        q(6) = -q(6);
        q(8) = q(8)+pi;
    end
    q(8) = phaseReference+atan2(sin(q(8)-phaseReference), ...
        cos(q(8)-phaseReference));
end


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



function result = runDampedGLSDC(q0,tm,y,sigma,opts,maxIter,costTol,stepTol)
    % Same GLSDC least-squares correction as Task 2.
    % Backtracking is a numerical safeguard: q_new = q + alpha*dq.
    % Start alpha at 1 and halve it until the nonlinear cost decreases.
    result = struct('q',q0,'converged',false,'status','Maximum iterations', ...
        'sol',[],'J',Inf,'backtracked',false);
    q = q0;
    previousJ = Inf;
    previousAlpha = 1;
    for iteration = 0:maxIter
        try
            [ypred,H,sol] = propagateModel(q,tm,opts);
        catch problem
            result.status = ['Integration failed: ' problem.message];
            return
        end
        e = y-ypred;
        J = 0.5*sum((e/sigma).^2);
        A = H/sigma;
        b = e/sigma;
        c = sqrt(sum(A.^2,1));
        if any(c==0) || any(~isfinite(c))
            result.status = 'Invalid sensitivity column';
            return
        end
        As = A./c;
        if rank(As)<8
            result.status = 'Measurement Jacobian lost rank';
            return
        end
        dq = (As\b)./c.';
        relativeStep = max(abs(dq)./max(1,abs(q)));
        relativeCost = abs(previousJ-J)/max(previousJ,eps);
        result.q = q;
        result.sol = sol;
        result.J = J;
        % A step-size guard prevents a tiny backtracking step from being
        % mistaken for meaningful cost convergence.
        costConverged = iteration>0 && relativeCost<costTol && ...
            previousAlpha>=1/16;
        if relativeStep<stepTol || costConverged
            result.converged = true;
            result.status = 'Converged';
            return
        end
        if iteration==maxIter
            return
        end
        accepted = false;
        alpha = 1;
        for backtrack = 0:20
            qtrial = q+alpha*dq;
            Jtrial = positionOnlyCost(qtrial,tm,y,sigma,opts);
            if isfinite(Jtrial) && Jtrial<=J+1e-11*max(1,J)
                accepted = true;
                break
            end
            alpha = alpha/2;
        end
        if ~accepted
            result.status = 'No decreasing GLSDC step found';
            return
        end
        result.backtracked = result.backtracked || alpha<1;
        previousJ = J;
        previousAlpha = alpha;
        q = qtrial;
    end
end

function J = positionOnlyCost(q,tm,y,sigma,opts)
    J = Inf;
    if any(~isfinite(q))
        return
    end
    try
        sol = ode45(@(t,s) stateOnlyODE(t,s,q(3:8)),[0 tm(end)],q(1:2),opts);
        if sol.x(end)<tm(end)-1e-9
            return
        end
        S = deval(sol,tm);
        residual = y-S(1,:).';
        if all(isfinite(residual))
            J = 0.5*sum((residual/sigma).^2);
        end
    catch
        % A failed trial step is rejected without changing the measurements.
    end
end

function ds = stateOnlyODE(t,s,p)
    x = s(1);
    v = s(2);
    ds = [v;-p(1)*v-p(2)*x-p(3)*x^3-p(4)*sin(p(5)*t+p(6))];
end

