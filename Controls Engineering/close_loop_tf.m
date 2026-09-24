% Plant
A = [0 1; 0 -10];  B = [0;1];  C = [1 0];  D = 0;

% Specs → controller
zeta = 0.5; wn = 85;
K = [wn^2, 2*zeta*wn - 10];   % -> [7225, 75]
Acl = A - B*K;
Nbar = -1/(C*(Acl\B));        % -> 7225

% Closed-loop (reference → y)
sys_cl = ss(Acl, B*Nbar, C, D);
step(sys_cl), grid on
S = stepinfo(sys_cl);       
eig(Acl)                     

