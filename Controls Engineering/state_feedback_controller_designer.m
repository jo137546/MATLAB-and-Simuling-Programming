% Plant
A = [0 1; 0 -10];
B = [0; 1];
C = [1 0];
D = 0;

% Controller from specs
zeta = 0.404; wn = 100;
K = [wn^2, 2*zeta*wn - 10];     % k - gain  
Acl = A - B*K;
Nbar = -1/(C*(Acl\B));          

% Closed-loop (reference -> y)
sys_cl = ss(Acl, B*Nbar, C, 0);
S = stepinfo(sys_cl);            % check Overshoot and SettlingTime
disp(K); disp(Nbar); disp(S);
figure; step(sys_cl), grid on
title('Closed-loop step (ref \rightarrow y)');