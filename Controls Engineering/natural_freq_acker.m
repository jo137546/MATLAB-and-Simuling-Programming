%% Pole placement with 'acker' for a 2x2 SISO system
% Requires Control System Toolbox.

% --- GIVEN SYSTEM ---
A = [0 1;
    -6 -5];
B = [0;
     1];

% --- SPECS (edit these) ---
zeta = 0.707;     % damping ratio
Tp   = 3.14;      % peak time (seconds)

% --- Compute desired natural and damped frequencies ---
wd = pi / Tp;                         % damped natural frequency
wn = wd / sqrt(1 - zeta^2);           % natural frequency

% Desired 2nd-order closed-loop polynomial: s^2 + 2*zeta*wn*s + wn^2
% Desired poles (complex conjugates):
p = [-zeta*wn + 1i*wn*sqrt(1 - zeta^2), ...
     -zeta*wn - 1i*wn*sqrt(1 - zeta^2)];

% --- Controllability check ---
Co = ctrb(A,B);
if rank(Co) < size(A,1)
    error('System is NOT controllable; pole placement impossible.');
end

% --- Compute state-feedback gain using ACKER ---
K = acker(A,B,p);

% --- Verify poles ---
eigs_cl = eig(A - B*K);

% --- Display results ---
fprintf('Specs:  zeta = %.4f,  Tp = %.4f s\n', zeta, Tp);
fprintf('wn = %.4f rad/s,  wd = %.4f rad/s\n', wn, wd);
fprintf('Desired poles: %s\n', mat2str(p,4));
fprintf('K = %s\n', mat2str(K,4));
fprintf('eig(A-B*K) = %s\n', mat2str(eigs_cl,4));

% --- OPTIONAL quick step check (u is the input, no prefilter) ---
% C = [1 0]; D = 0;
% sys_cl = ss(A - B*K, B, C, D);
% figure; step(sys_cl); grid on; title('Closed-loop (input u to output y)');