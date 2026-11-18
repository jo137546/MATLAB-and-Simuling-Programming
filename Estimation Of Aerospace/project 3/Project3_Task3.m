clear; clc;

base = baselineTask2();
epochs = base.epochs;
qEstimate = base.qEstimate;
P0 = base.P0;
[Pz,epochErrors,T] = propagateErrorCovariance(qEstimate,P0, ...
    base.solEstimate,epochs,base.qTrue,base.trueEpoch);
sigmaAtEpochs = zeros(8,4);
jointDistanceSquared = zeros(4,1);
for k = 1:4
    sigmaAtEpochs(:,k) = sqrt(diag(Pz(:,:,k)));
    e = epochErrors(1:2,k);
    jointDistanceSquared(k) = e.'*(Pz(1:2,1:2,k)\e);
end
absErrorOverSigma = abs(epochErrors)./sigmaAtEpochs;
comparisonTable = table(kron(epochs(:),ones(8,1)),repmat(base.names,4,1), ...
    epochErrors(:),sigmaAtEpochs(:),absErrorOverSigma(:), ...
    absErrorOverSigma(:)<=1,absErrorOverSigma(:)<=2,absErrorOverSigma(:)<=3, ...
    'VariableNames',{'Time','Quantity','ActualError','OneSigma', ...
    'AbsErrorOverSigma','Within1Sigma','Within2Sigma','Within3Sigma'});
fprintf('\nTask 3: covariance of initial states and parameters\n');
disp(P0);
fprintf('Task 3: all eight errors and marginal standard deviations\n');
disp(comparisonTable);
jointEllipseTable = table(epochs(:),jointDistanceSquared, ...
    jointDistanceSquared<=1,jointDistanceSquared<=4,jointDistanceSquared<=9, ...
    'VariableNames',{'Time','SquaredMahalanobisDistance','Inside1Sigma', ...
    'Inside2Sigma','Inside3Sigma'});
disp(jointEllipseTable);
for k = 1:4
    fprintf('\nState covariance at t = %g\n',epochs(k));
    disp(Pz(1:2,1:2,k));
end

figure('Name','Task 3 - Position and velocity error ellipses', ...
    'NumberTitle','off','Visible','on','Color','w','Position',[80 80 1150 850]);
tiledlayout(2,2,'TileSpacing','compact');
for k = 1:4
    nexttile; hold on;
    h = drawEllipses(Pz(1:2,1:2,k),[0;0],[1 2 3],'-');
    hp = plot(epochErrors(1,k),epochErrors(2,k),'kx', ...
        'MarkerSize',9,'LineWidth',1.8);
    xline(0,'k:','HandleVisibility','off');
    yline(0,'k:','HandleVisibility','off');
    grid on; axis equal;
    xlabel('Position error, estimated - true');
    ylabel('Velocity error, estimated - true');
    title(sprintf('t = %g',epochs(k)));
    legend([h(:);hp],{'1 sigma','2 sigma','3 sigma','Actual error'}, ...
        'Location','best');
end
fprintf(['\nThese are ellipses in ERROR coordinates, centered at zero.\n' ...
    'For a joint 2-D Gaussian, radii 1, 2, 3 contain approximately ' ...
    '39.35%%, 86.47%%, 98.89%% probability.\n' ...
    'A marginal error/sigma comparison and joint ellipse membership ' ...
    'are different checks.\n']);
% Structural checks on the full 8-by-8 propagation.
assert(norm(T(:,:,1)-eye(8),'fro')<1e-10,'T(0) should be identity.');
assert(norm(Pz(:,:,1)-P0,'fro')<1e-10*max(1,norm(P0,'fro')));
for k = 1:4
    assert(norm(Pz(3:8,3:8,k)-P0(3:8,3:8),'fro')< ...
        1e-10*max(1,norm(P0,'fro')),'Constant-parameter covariance changed.');
end
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


