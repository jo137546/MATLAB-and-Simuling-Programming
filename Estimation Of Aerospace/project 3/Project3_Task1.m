clc; clear; close all
rng(6414,'twister');           % Reproducible measurement-noise realization

%% 1. Given values and initial conditions
p = [0.05; 4; 0.2; -0.5; 10; pi/2]; % IMPORTANT: p4 is negative
s0 = [2; 0];
sigma = 0.1;
tm = (1:3000)'/10;            % Measurements start at 0.1, NOT at zero
t = (0:0.02:300)';            % Dense output grid for smooth motion plots
Phi0 = eye(2);
Psi0 = zeros(2,6);
z0 = [s0; Phi0(:); Psi0(:)];   % 2 + 4 + 12 = 18 entries
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',0.03);

%% 2. Integrate motion and both sensitivity matrices together
sol = ode45(@(t,z) augmentedODE(t,z,p),[0 300],z0,opts);
Z = deval(sol,t).';           % Each row corresponds to one output time
Zm = deval(sol,tm).';
S = Z(:,1:2);
Phi = reshape(Z(:,3:6).',2,2,[]);
Psi = reshape(Z(:,7:18).',2,6,[]);
Phi_m = reshape(Zm(:,3:6).',2,2,[]);
Psi_m = reshape(Zm(:,7:18).',2,6,[]);
xtrue_m = Zm(:,1);

%% 3. Generate position measurements (no process noise)
noise = sigma*randn(numel(tm),1); % randn has variance 1; scale by sigma
ymeas = xtrue_m + noise;
fprintf('Project 3, Task 1\n');
fprintf('p = [%.6g %.6g %.6g %.6g %.6g %.12g]\n',p);
fprintf('Number of measurements: %d\n',numel(tm));
fprintf('Measurement interval: %.1f to %.1f\n',tm(1),tm(end));
fprintf('Final state at t=300: x=%.12g, xdot=%.12g\n',S(end,:));
fprintf('Noise mean: %.8g (target 0)\n',mean(noise));
fprintf('Noise sample variance: %.8g (target %.8g)\n',var(noise,0),sigma^2);
fprintf('Initial acceleration: %.8g\n',stateODE(0,s0,p)'*[0;1]);

%% 4. Validate ALL columns of Phi and Psi using centered differences
% These are independent integrations of ONLY the two physical state ODEs.
% h = eta*max(1,abs(q)). Check two perturbation sizes to assess sensitivity
% to finite-difference step. Relative errors use a whole-history norm;
% they do not divide point by point by entries that may cross zero.
etaValues = [1e-4 1e-5];
q = [s0; p];
labels = {'x0'; 'xdot0'; 'p1'; 'p2'; 'p3'; 'p4'; 'p5'; 'p6'};
optsFD = odeset('RelTol',1e-11,'AbsTol',1e-13,'MaxStep',0.03);
relativeErrors = zeros(8,numel(etaValues));
absoluteErrors = zeros(8,numel(etaValues));
stepSizes = zeros(8,numel(etaValues));
for e = 1:numel(etaValues)
    FD = zeros(2,8,numel(tm));
    for k = 1:8
        h = etaValues(e)*max(1,abs(q(k)));
        stepSizes(k,e) = h;
        qp = q; qm = q;
        qp(k) = qp(k)+h;
        qm(k) = qm(k)-h;
        solPlus = ode45(@(t,s) stateODE(t,s,qp(3:8)), ...
                       [0 300],qp(1:2),optsFD);
        solMinus = ode45(@(t,s) stateODE(t,s,qm(3:8)), ...
                        [0 300],qm(1:2),optsFD);
        fd = (deval(solPlus,tm)-deval(solMinus,tm))/(2*h);
        FD(:,k,:) = reshape(fd,2,1,[]);
        if k <= 2
            reference = reshape(Phi_m(:,k,:),2,[]);
        else
            reference = reshape(Psi_m(:,k-2,:),2,[]);
        end
        error = fd-reference;
        relativeErrors(k,e) = norm(error,'fro')/max(norm(reference,'fro'),eps);
        absoluteErrors(k,e) = max(abs(error),[],'all');
    end
    fprintf('\nCentered-difference validation: eta = %.0e\n',etaValues(e));
    T = table(labels,stepSizes(:,e),relativeErrors(:,e),absoluteErrors(:,e), ...
        'VariableNames',{'Quantity','Perturbation','RelativeHistoryError','MaxAbsoluteError'});
    disp(T);
end
% Retain the smaller-perturbation histories for the validation plots.
PhiFD = FD(:,1:2,:);
PsiFD = FD(:,3:8,:);
% Ordinary signed error for every matrix entry: ODE minus finite difference.
% Keep individual entries instead of combining them into a matrix norm.
errPhi = Phi_m-PhiFD;          % 2 x 2 x number of measurement times
errPsi = Psi_m-PsiFD;          % 2 x 6 x number of measurement times

% Additional exact check: trace(F)=-p1, hence det(Phi)=exp(-p1*t).
detPhi = squeeze(Phi(1,1,:).*Phi(2,2,:)-Phi(1,2,:).*Phi(2,1,:));
detExact = exp(-p(1)*t);
detRelativeError = max(abs(detPhi-detExact)./detExact);
fprintf('Max relative determinant identity error: %.8g\n',detRelativeError);
fprintf('Phi(0) identity error: %.8g\n',norm(Phi(:,:,1)-eye(2),'fro'));
fprintf('Psi(0) zero error: %.8g\n',norm(Psi(:,:,1),'fro'));
assert(all(isfinite(Z),'all'),'Nonfinite trajectory.');
assert(all(relativeErrors(:,end)<1e-3),'Inspect finite-difference validation.');
assert(detRelativeError<1e-4,'Inspect state-transition determinant check.');

%% 5. Required plots
f = figure('Name','True states','NumberTitle','off','Visible','on', ...
    'Color','w','Position',[100 100 1100 700]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile; plot(t,S(:,1),'LineWidth',1); grid on;
xlabel('Time, t'); ylabel('Position, x'); title('True position');
nexttile; plot(t,S(:,2),'LineWidth',1); grid on;
xlabel('Time, t'); ylabel('Velocity, dx/dt'); title('True velocity');

f = figure('Name','Phase portrait','NumberTitle','off','Visible','on', ...
    'Color','w','Position',[100 100 1000 500]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile; plot(S(:,1),S(:,2),'LineWidth',0.8); grid on;
xlabel('Position, x'); ylabel('Velocity, dx/dt'); title('Phase portrait: 0 to 300');
nexttile; tail = t>=290;
plot(S(tail,1),S(tail,2),'LineWidth',1); grid on;
xlabel('Position, x'); ylabel('Velocity, dx/dt'); title('Late response: 290 to 300');

f = figure('Name','Simulated measurements','NumberTitle','off','Visible','on', ...
    'Color','w','Position',[100 100 1100 700]);
tiledlayout(2,1,'TileSpacing','compact');
nexttile;
plot(tm,ymeas,'.','Color',[0.7 0.7 0.7],'MarkerSize',4); hold on;
plot(t,S(:,1),'b','LineWidth',1); grid on;
xlabel('Time, t'); ylabel('Position'); title('Measurements and true motion');
legend('Measurements','True position');
nexttile;
plot(tm,ymeas,'.','Color',[0.5 0.5 0.5],'MarkerSize',7); hold on;
plot(t,S(:,1),'b','LineWidth',1.2); grid on; xlim([290 300]);
xlabel('Time, t'); ylabel('Position'); title('Late response with measurement noise');
legend('Measurements','True position');

%% 6. Ordinary sensitivity error plots (linear axes, zero is exact agreement)
f = figure('Name','Phi errors','NumberTitle','off','Visible','on', ...
    'Color','w','Position',[100 100 1100 700]);
tiledlayout(2,2,'TileSpacing','compact');
for row = 1:2
    for col = 1:2
        nexttile;
        plot(tm,squeeze(errPhi(row,col,:)),'b','LineWidth',1);
        yline(0,'k:'); grid on;
        xlabel('Time, t'); ylabel('ODE - finite difference');
        title(sprintf('Error in \\Phi_{%d%d}',row,col));
    end
end
sgtitle('State-transition sensitivity errors');

% Separate position and velocity sensitivity errors for readable scales.
stateNames = {'Position','Velocity'};
for row = 1:2
    f = figure('Name',['Psi errors - ' stateNames{row}], ...
        'NumberTitle','off','Visible','on','Color','w', ...
        'Position',[100 100 1200 750]);
    tiledlayout(2,3,'TileSpacing','compact');
    for col = 1:6
        nexttile;
        plot(tm,squeeze(errPsi(row,col,:)),'b','LineWidth',1);
        yline(0,'k:'); grid on;
        xlabel('Time, t'); ylabel('ODE - finite difference');
        title(sprintf('Error in \\Psi_{%d%d} (p_%d)',row,col,col));
    end
    sgtitle([stateNames{row} ' parameter-sensitivity errors']);
end

drawnow;
fprintf('\nPlots are open in MATLAB figures. Numerical data remain in the Workspace.\n');

%% Local functions
function ds = stateODE(t,s,p)
    x = s(1);
    v = s(2);
    ds = [v; -p(1)*v-p(2)*x-p(3)*x^3-p(4)*sin(p(5)*t+p(6))];
end

function dz = augmentedODE(t,z,p)
    s = z(1:2);
    x = s(1);
    v = s(2);
    Phi = reshape(z(3:6),2,2);
    Psi = reshape(z(7:18),2,6);
    theta = p(5)*t+p(6);
    % F = partial f / partial s, evaluated along the current trajectory.
    F = [0, 1; -p(2)-3*p(3)*x^2, -p(1)];
    % G = partial f / partial p (holding the current state fixed).
    G = [zeros(1,6); ...
        -v, -x, -x^3, -sin(theta), -p(4)*t*cos(theta), -p(4)*cos(theta)];
    dPhi = F*Phi;
    dPsi = F*Psi+G;
    dz = [stateODE(t,s,p); dPhi(:); dPsi(:)];
end

